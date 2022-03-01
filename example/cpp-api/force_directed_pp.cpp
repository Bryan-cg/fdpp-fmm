#include "force_directed_pp.hpp"

using namespace FDPP;
using namespace FDPP::FORCE;

#define EARTH_R 6378137
#define d2r (M_PI / 180.0)             // degrees to radians
#define r2d (180.0 / M_PI)             // radians to degrees
#define a2m ((M_PI / 180.0) * EARTH_R) // angular units to meters
#define Haversine bg::strategy::distance::haversine<double>(EARTH_R)

ElectricForce::ElectricForce(double length_edge_arg, double distance_arg, double c_arg) : length_edge(length_edge_arg), distance(distance_arg),
                                                                                          c(c_arg){};

double ElectricForce::calculate(double cos_theta) const
{
    double force;
    if (cos_theta < 0)
    {
        cos_theta = cos_theta / 2;
    }
    if (distance <= 0.000)
    {
        force = 0; 
    }
    else
    {
        force = 1 * sqrt(length_edge) * cos_theta / distance;
    }
    return force;
};

double ElectricForce::calculate_cos_theta(const Point &p1_edge, const Point &p2_edge, const Point &p1_trace, const Point &p2_trace) const
{
    // Convert points to radians for calculations formula
    double az_edge = bg::formula::spherical_azimuth(p1_edge.get<0>() * d2r, p1_edge.get<1>() * d2r, p2_edge.get<0>() * d2r, p2_edge.get<1>() * d2r);
    double az_trace = bg::formula::spherical_azimuth(p1_trace.get<0>() * d2r, p1_trace.get<1>() * d2r, p2_trace.get<0>() * d2r, p2_trace.get<1>() * d2r);
    return cos(az_trace - az_edge);
};

SpringForceAttr::SpringForceAttr(double distance_arg, double c1_arg, double c2_arg) : distance(distance_arg), c1(c1_arg), c2(c2_arg){};

double SpringForceAttr::calculate() const
{
    if (distance <= 0.000)
    {
        return 0.0;
    }
    return c1 * log(distance / c2);
};

SpringForceRep::SpringForceRep(double distance_arg, double c3_arg) : distance(distance_arg), c3(c3_arg){};

double SpringForceRep::calculate() const
{
    if (distance <= 0.000)
    {
        return 0.0;
    }
    return c3 / (distance * distance * distance);
};

ForceDirectedPP::ForceDirectedPP(char *shapeFile, char *tracesFile, char *ubodtFile, int iterations_arg)
{
    std::string shapeFile_str(shapeFile);
    std::string tracesFile_str(tracesFile);

    this->_iterations_fdpp = iterations_arg;
    this->network = new NETWORK::Network(shapeFile_str, "fid", "u", "v");
    this->networkGraph = new NETWORK::NetworkGraph(*(this->network));
    this->gps_config = new GPSConfig(tracesFile);
    this->result_config = new ResultConfig();
    result_config->file = "output_fdpp_5.csv";

    SPDLOG_INFO("Iterations FDPP {}", _iterations_fdpp);
    SPDLOG_INFO("Network configuration FDPP - FMM");
    SPDLOG_INFO("Network file {}", shapeFile_str);
    SPDLOG_INFO("Network node count {}", network->get_node_count());
    SPDLOG_INFO("Network edge count {}", network->get_edge_count());
    SPDLOG_INFO("NetworkGraph number of vertices {}", networkGraph->get_num_vertices());
    SPDLOG_INFO("Network initialized");
    gps_config->print();
    result_config->print();
    SPDLOG_INFO("Initializing FMM");
    fmmw = std::unique_ptr<fmm_wrap>(new fmm_wrap(ubodtFile, *network, *networkGraph));
}

ForceDirectedPP::~ForceDirectedPP()
{
    delete network;
    delete networkGraph;
    delete gps_config;
    delete result_config;
}

LineString ForceDirectedPP::point_to_lineString(const Point &p)
{
    LineString temp;
    temp.add_point(p);
    return temp;
}

double ForceDirectedPP::interpolated_distance_lineString(const LineString &ls)
{
    int N = ls.get_num_points();
    Point p1 = ls.get_point(0);
    Point p2 = ls.get_point(N - 1);
    return haversine_distance_m(p1, p2);
}

double ForceDirectedPP::haversine_distance_m(const Point &p1, const Point &p2)
{
    // Make degree points for distance calculations
    PointDeg _p1(p1.get<0>(), p1.get<1>());
    PointDeg _p2(p2.get<0>(), p2.get<1>());
    return bg::distance(_p1, _p2, Haversine);
}

void ForceDirectedPP::match()
{
    if (!result_config->validate())
    {
        SPDLOG_CRITICAL("result_config invalid");
        return;
    }
    FMM::IO::GPSReader reader(*gps_config);
    FMM::IO::CSVMatchResultWriter writer(result_config->file, result_config->output_config);
    int buffer_trajectories_size = 100000;
    while (reader.has_next_trajectory())
    {
        SPDLOG_INFO("Force directed displacement of buffered GPS points");
        std::vector<Trajectory> trajectories = reader.read_next_N_trajectories(buffer_trajectories_size);
        int trajectories_fetched = trajectories.size();
        for (int i = 0; i < trajectories_fetched; ++i)
        {
            Trajectory &trajectory = trajectories[i];
            force_directed_displacement(trajectory);
        }
        SPDLOG_INFO("FMM of displaced GPS points");
        std::vector<MM::MatchResult> match_results = fmmw->match(trajectories);
        SPDLOG_INFO("Flushing matched trajectories");
        for (int i = 0; i < trajectories_fetched; i++)
        {
            Trajectory &trajectory = trajectories[i];
            MatchResult &match_result = match_results[i];
            writer.write_result(trajectory, match_result);
        }
    }
    SPDLOG_INFO("FDPP & FMM completed");
}

void ForceDirectedPP::force_directed_displacement(Trajectory &trajectory)
{
    LineString ls = trajectory.geom;
    LineStringDeg ls_d = cart_to_degr_linestring(ls);
    // if (bg::length(ls_d) > 300 && ls.get_num_points() >= 30)
    // {
        // number of iterations
        for (int i = 0; i < 20; i++)
        {
            // Problem, only one direction for each edge in stead of two?
            for (int i = 1; i < ls.get_num_points() - 1; i++)
            {
                Point p = ls.get_point(i);
                LineString ls_point = point_to_lineString(p);

                // get edges within radius
                FMM::MM::Traj_Candidates candidates = network->search_tr_cs(ls_point, 0.001); // todo point

                // displacement point
                Point p_res = calculate_net_force(i, candidates, ls);
                ls.get_geometry().at(i).set<0>(p_res.get<0>());
                ls.get_geometry().at(i).set<1>(p_res.get<1>());
            }
        // }
    }
    trajectory.geom = ls;
    std::cout << ls << std::endl;
}

Point ForceDirectedPP::calculate_net_force(const int point_i, const FMM::MM::Traj_Candidates &candidates, const LineString &trace_ls)
{
    const Point p = trace_ls.get_point(point_i);
    const Point p_prev = trace_ls.get_point(point_i - 1);
    const Point p_next = trace_ls.get_point(point_i + 1);

    // calculate total attractive spring force
    std::vector<double> p_cart = to_cart_cord(p);
    std::vector<double> p1_delta = calc_attr_spring_force_displacement(p, p_prev);//different direction?
    std::vector<double> p2_delta = calc_attr_spring_force_displacement(p, p_next);
    double dx_s_attr = p1_delta[0]; //+ p2_delta[0];
    double dy_s_attr = p1_delta[1]; //+ p2_delta[1];
    double dz_s_attr = p1_delta[2]; //+ p2_delta[2];

    // calculate total repelling spring force displacement
    double dx_s_rep = 0.0;
    double dy_s_rep = 0.0;
    double dz_s_rep = 0.0;
    /*
    for (int i = 0; i < trace_ls.get_num_points(); i++)
    {
        // Not an adjacent vertex
        if (i != point_i - 1 && i != point_i && i != point_i + 1)
        {
            const Point non_adj_point = trace_ls.get_point(i);
            std::vector<double> non_adj_delta = calc_rep_spring_force_displacement(p, non_adj_point);
            dx_s_rep += non_adj_delta[0];
            dy_s_rep += non_adj_delta[1];
            dz_s_rep += non_adj_delta[2];
        }
    }
    */

    // calculate total electric force
    double fe_total = 0.0;
    std::vector<FMM::MM::Candidate> candidates_point_i;
    if (candidates.size() > 0)
    {
        candidates_point_i = candidates[0];
        std::vector<Point> force_directions;
        for (FMM::MM::Candidate j : candidates_point_i)
        {
            LineString edge_lst = j.edge->geom;
            LineStringDeg lsdg = cart_to_degr_linestring(edge_lst);
            Point p1_edge = edge_lst.get_point(0);
            Point p2_edge = edge_lst.get_point(edge_lst.get_num_points() - 1);
            Point p_dir;
            double d;
            double l = bg::length(lsdg); // haversine_distance_m(p1_edge, p2_edge);
            calculate_closest_point(edge_lst, p, &d, &p_dir);
            ElectricForce ef = {l, d};
            double cos_theta = ef.calculate_cos_theta(p1_edge, p2_edge, p_prev, p_next);
            double ef_i = ef.calculate(cos_theta);
            fe_total += ef_i;
            force_directions.push_back(p_dir);
        }

        double dx_e = 0.0;
        double dy_e = 0.0;
        double dz_e = 0.0;
        for (Point p_dir : force_directions)
        {
            std::vector<double> delta = calc_electric_force_displacement(p, p_dir, fe_total);
            dx_e += delta[0];
            dy_e += delta[1];
            dz_e += delta[2];
        }

        // Calculate replaced point
        std::vector<double> delta_res;
        delta_res.push_back(p_cart[0] + dx_s_attr - dx_s_rep + dx_e);
        delta_res.push_back(p_cart[1] + dy_s_attr - dy_s_rep + dy_e);
        delta_res.push_back(p_cart[2] + dz_s_attr - dz_s_rep + dz_e);
        Point p_res = from_cart_cord(delta_res);
        return p_res;
    }
    else
    {
        std::vector<double> delta_res;
        delta_res.push_back(p_cart[0] + dx_s_attr - dx_s_rep);
        delta_res.push_back(p_cart[1] + dy_s_attr - dy_s_rep);
        delta_res.push_back(p_cart[2] + dz_s_attr - dz_s_rep);
        Point p_res = from_cart_cord(delta_res);
        return p_res;
    }
}

std::vector<double> ForceDirectedPP::calc_attr_spring_force_displacement(const Point &p1, const Point &p2)
{
    std::vector<double> p1_cart = to_cart_cord(p1);
    std::vector<double> p2_cart = to_cart_cord(p2);

    double distance = haversine_distance_m(p1, p2);
    // Duplicated point
    if (distance <= 0.000)
    {
        std::vector<double> p_res_cart;
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        return p_res_cart;
    }
    SpringForceAttr sf_prev = {distance};
    double fs_total = sf_prev.calculate();
    double dx = abs(p1_cart[0] - p2_cart[0]);
    double dy = abs(p1_cart[1] - p2_cart[1]);
    double dz = abs(p1_cart[2] - p2_cart[2]);
    double fx = fs_total * (dx / distance);
    double fy = fs_total * (dy / distance);
    double fz = fs_total * (dz / distance);
    if (p2_cart[0] < p1_cart[0])
        fx = -fx;
    if (p2_cart[1] < p1_cart[1])
        fy = -fy;
    if (p2_cart[2] < p1_cart[2])
        fz = -fz;

    double force_dist = sqrt(fx * fx + fy * fy + fz * fz);
    double max = 80; // 80m
    double limitedDist = std::min(force_dist, max);
    double dx_x = 0.2 * (fx / force_dist * limitedDist);
    double dy_y = 0.2 * (fy / force_dist * limitedDist);
    double dz_z = 0.2 * (fz / force_dist * limitedDist);
    std::vector<double> p_res_cart;
    p_res_cart.push_back(dx_x);
    p_res_cart.push_back(dy_y);
    p_res_cart.push_back(dz_z);
    return p_res_cart;
}

std::vector<double> ForceDirectedPP::calc_rep_spring_force_displacement(const Point &p1, const Point &p2)
{
    std::vector<double> p1_cart = to_cart_cord(p1);
    std::vector<double> p2_cart = to_cart_cord(p2);

    double distance = haversine_distance_m(p1, p2);
    // Duplicated point
    if (distance <= 0.000)
    {
        std::vector<double> p_res_cart;
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        return p_res_cart;
    }
    SpringForceRep sf_prev = {distance};
    double fs_total = sf_prev.calculate();
    double dx = abs(p1_cart[0] - p2_cart[0]);
    double dy = abs(p1_cart[1] - p2_cart[1]);
    double dz = abs(p1_cart[2] - p2_cart[2]);
    double fx = fs_total * (dx / distance);
    double fy = fs_total * (dy / distance);
    double fz = fs_total * (dz / distance);
    if (p2_cart[0] < p1_cart[0])
        fx = -fx;
    if (p2_cart[1] < p1_cart[1])
        fy = -fy;
    if (p2_cart[2] < p1_cart[2])
        fz = -fz;

    double force_dist = sqrt(fx * fx + fy * fy + fz * fz);
    double max = 80; // 80m
    double limitedDist = std::min(force_dist, max);
    double dx_x = 0.0001 * (fx / force_dist * limitedDist);
    double dy_y = 0.0001 * (fy / force_dist * limitedDist);
    double dz_z = 0.0001 * (fz / force_dist * limitedDist);
    std::vector<double> p_res_cart;
    p_res_cart.push_back(dx_x);
    p_res_cart.push_back(dy_y);
    p_res_cart.push_back(dz_z);
    return p_res_cart;
}

std::vector<double> ForceDirectedPP::calc_electric_force_displacement(const Point &p1, const Point &p2, const double total_force)
{
    std::vector<double> p1_cart = to_cart_cord(p1);
    std::vector<double> p2_cart = to_cart_cord(p2);
    double distance = haversine_distance_m(p1, p2);
    double dx = abs(p1_cart[0] - p2_cart[0]);
    double dy = abs(p1_cart[1] - p2_cart[1]);
    double dz = abs(p1_cart[2] - p2_cart[2]);
    double fx = total_force * (dx / distance);
    double fy = total_force * (dy / distance);
    double fz = total_force * (dz / distance);
    if (p2_cart[0] < p1_cart[0])
        fx = -fx;
    if (p2_cart[1] < p1_cart[1])
        fy = -fy;
    if (p2_cart[2] < p1_cart[2])
        fz = -fz;

    double force_dist = sqrt(fx * fx + fy * fy + fz * fz);
    double max = 80; // 80m
    double limitedDist = std::min(force_dist, max);
    double dx_x = 0.01 * (fx / force_dist * limitedDist);
    double dy_y = 0.01 * (fy / force_dist * limitedDist);
    double dz_z = 0.01 * (fz / force_dist * limitedDist);
    std::vector<double> p_res_cart;
    p_res_cart.push_back(dx_x);
    p_res_cart.push_back(dy_y);
    p_res_cart.push_back(dz_z);
    return p_res_cart;
}

// x, y, z
std::vector<double> ForceDirectedPP::to_cart_cord(const Point &p)
{
    std::vector<double> result;
    result.push_back(EARTH_R * cos(p.get<1>() * d2r) * cos(p.get<0>() * d2r));
    result.push_back(EARTH_R * cos(p.get<1>() * d2r) * sin(p.get<0>() * d2r));
    result.push_back(EARTH_R * sin(p.get<1>() * d2r));
    return result;
}

Point ForceDirectedPP::from_cart_cord(const std::vector<double> &coord)
{
    double lat = r2d * (asin(coord[2] / EARTH_R));
    double lon = r2d * (atan2(coord[1], coord[0]));
    Point p;
    p.set<0>(lon);
    p.set<1>(lat);
    return p;
}

void ForceDirectedPP::calculate_closest_point(const LineString &ls_edge, const Point &p, double *dist, Point *c_p)
{
    double px = p.get<0>();
    double py = p.get<1>();
    double offset;
    double closest_x, closest_y;
    ALGORITHM::linear_referencing(px, py, ls_edge, dist, &offset, &closest_x, &closest_y);
    *c_p = Point(closest_x, closest_y);
    *dist = haversine_distance_m(p, *c_p);
}

LineStringDeg ForceDirectedPP::cart_to_degr_linestring(const LineString &ls)
{
    LineStringDeg lsdg;
    for (int i = 0; i < ls.get_num_points(); i++)
    {
        Point p = ls.get_point(i);
        PointDeg pdg(p.get<0>(), p.get<1>());
        bg::append(lsdg, pdg);
    }
    return lsdg;
}
