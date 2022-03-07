#include "force_directed_pp.hpp"

using namespace FDPP;
using namespace FDPP::FORCE;

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
        force = c * sqrt(length_edge) * cos_theta;
    }
    else
    {
        force = c * sqrt(length_edge) * cos_theta / distance;
    }
    return force;
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
    result_config->file = "output_fdpp_10.csv";

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
        std::vector<Trajectory> trajectories_init = reader.read_next_N_trajectories(buffer_trajectories_size);
        std::vector<Trajectory> trajectories;
        for (auto &&t : trajectories_init)
        {
            Trajectory tn(t);
            tn.geom = UTIL::lower_sample_freq(t.geom, 15);
            trajectories.push_back(tn);
        }
        std::cout << trajectories_init[2].geom << std::endl;
        std::cout << trajectories[2].geom << std::endl;
        std::vector<Trajectory> traces = trajectories;
        int trajectories_fetched = trajectories.size();
        for (int i = 0; i < trajectories_fetched; ++i)
        {
            Trajectory &trajectory = trajectories[i];
            force_directed_displacement(trajectory);
        }
        SPDLOG_INFO("FMM of displaced GPS points");
        std::vector<MM::MatchResult> mr_pp = fmmw->match(trajectories);
        SPDLOG_INFO("FMM of original trajectories");
        std::vector<MM::MatchResult> mr_no_pp = fmmw->match(traces);
        sort(mr_pp.begin(), mr_pp.end(), [](const MatchResult &r1, const MatchResult &r2)
             { return r1.id < r2.id; });
        sort(mr_no_pp.begin(), mr_no_pp.end(), [](const MatchResult &r1, const MatchResult &r2)
             { return r1.id < r2.id; });

        double count_improved = 0.0;
        std::vector<MM::MatchResult> final_results = combine_fmm_fdpp_output(mr_pp, mr_no_pp, traces, &count_improved);

        SPDLOG_INFO("Calculating accuracy & flushing matches to output file");

        double li_total = 0.0;
        double avgE_total = 0.0;
        double frechet_total = 0.0;
        double hausdorff_total = 0.0;
        double count_unmatched = 0.0;
        for (int i = 0; i < trajectories_fetched; i++)
        {
            Trajectory &trajectory = trajectories[i];
            MatchResult &match_result = final_results[i];
            writer.write_result(trajectory, match_result);
            double li, avgE, frechet, hausdorff;
            GEOM::calc_accuracy(
                GEOM::cart_to_degr_linestring(traces[i].geom),
                GEOM::cart_to_degr_linestring(match_result.mgeom),
                &li,
                &avgE,
                &frechet,
                &hausdorff);
            if (li == INT_MAX)
            {
                count_unmatched += 1.0;
                continue;
            }
            li_total += li;
            avgE_total += avgE;
            frechet_total += frechet;
            hausdorff_total += hausdorff;
        }

        SPDLOG_INFO("Accuracy: li {}, avg error {}, frechet {}, hausdorff {}",
                    li_total / trajectories_fetched,
                    avgE_total / trajectories_fetched,
                    frechet_total / trajectories_fetched,
                    hausdorff_total / trajectories_fetched);
        double improved_per = (count_improved / trajectories_fetched) * 100;
        SPDLOG_INFO("improved {}, total {}, percentage {}", count_improved, trajectories_fetched, improved_per);
        SPDLOG_INFO("unmatched traces {}", count_unmatched);
    }
    SPDLOG_INFO("FDPP & FMM completed");
}

std::vector<MM::MatchResult> ForceDirectedPP::combine_fmm_fdpp_output(
    const std::vector<MM::MatchResult> &mr_pp,
    const std::vector<MM::MatchResult> &mr_no_pp,
    const std::vector<Trajectory> &traces,
    double *count_improved)
{
    SPDLOG_INFO("Calculating most accurate matches");
    std::vector<MM::MatchResult> final_results;
    int trajectories_fetched = traces.size();
    for (int i = 0; i < trajectories_fetched; i++)
    {
        double li, avgE, frechet, hausdorff;
        double li_fmm, avgE_fmm, frechet_fmm, hausdorff_fmm;
        GEOM::calc_accuracy(
            GEOM::cart_to_degr_linestring(traces[i].geom),
            GEOM::cart_to_degr_linestring(mr_no_pp[i].mgeom),
            &li_fmm,
            &avgE_fmm,
            &frechet_fmm,
            &hausdorff_fmm);

        GEOM::calc_accuracy(
            GEOM::cart_to_degr_linestring(traces[i].geom),
            GEOM::cart_to_degr_linestring(mr_pp[i].mgeom),
            &li,
            &avgE,
            &frechet,
            &hausdorff);

        if (frechet < frechet_fmm && hausdorff < hausdorff_fmm)
        {
            *count_improved += 1.0;
            final_results.push_back(mr_pp[i]);
            // std::cout << traces[i].geom << std::endl;
            // std::cout << mr_pp[i].mgeom << std::endl;
            // std::cout << mr_no_pp[i].mgeom << std::endl;
        }
        else
        {
            final_results.push_back(mr_no_pp[i]);
        }
    }
    return final_results;
}

void ForceDirectedPP::force_directed_displacement(Trajectory &trajectory)
{
    FDPP::IO::CSVIterationWriter iter_writer("iterations.csv");
    LineString ls = trajectory.geom;
    LineStringDeg ls_d = GEOM::cart_to_degr_linestring(ls);

    // number of iterations
    for (int i = 0; i < _iterations_fdpp; i++)
    {
        // Problem, only one direction for each edge in stead of two?
        for (int j = 1; j < ls.get_num_points() - 1; j++)
        {
            Point p = ls.get_point(j);
            LineString ls_point = GEOM::point_to_lineString(p);

            // get edges within radius
            FMM::MM::Traj_Candidates candidates = network->search_tr_cs(ls_point, 0.001); // todo point

            // displacement point
            Point p_res = calculate_net_force(j, candidates, ls);
            ls.get_geometry().at(j).set<0>(p_res.get<0>());
            ls.get_geometry().at(j).set<1>(p_res.get<1>());
        }
        // iter_writer.write_result(i + 1, ls);
    }

    trajectory.geom = ls;
    // std::cout << ls << std::endl;
}

Point ForceDirectedPP::calculate_net_force(const int point_i, const FMM::MM::Traj_Candidates &candidates, const LineString &trace_ls)
{
    const Point p = trace_ls.get_point(point_i);
    const Point p_prev = trace_ls.get_point(point_i - 1);
    const Point p_next = trace_ls.get_point(point_i + 1);

    // calculate total attractive spring force
    std::vector<double> p_cart = GEOM::to_cart_cord(p);
    std::vector<double> p1_delta = calc_attr_spring_force_displacement(p, p_prev); // different direction?
    std::vector<double> p2_delta = calc_attr_spring_force_displacement(p, p_next);
    double dx_s_attr = p1_delta[0] + p2_delta[0];
    double dy_s_attr = p1_delta[1] + p2_delta[1];
    double dz_s_attr = p1_delta[2] + p2_delta[2];

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
            LineStringDeg lsdg = GEOM::cart_to_degr_linestring(edge_lst);
            Point p1_edge = edge_lst.get_point(0);
            Point p2_edge = edge_lst.get_point(edge_lst.get_num_points() - 1);
            Point p_dir;
            double d;
            double l = bg::length(lsdg); // haversine_distance_m(p1_edge, p2_edge);
            GEOM::calculate_closest_point(edge_lst, p, &d, &p_dir);
            ElectricForce ef = {l, d};
            double cos_theta = GEOM::calculate_cos_theta(p1_edge, p2_edge, p_prev, p_next);
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
        delta_res.push_back(p_cart[0] + dx_s_attr + dx_e);
        delta_res.push_back(p_cart[1] + dy_s_attr + dy_e);
        delta_res.push_back(p_cart[2] + dz_s_attr + dz_e);
        Point p_res = GEOM::from_cart_cord(delta_res);
        return p_res;
    }
    else
    {
        std::vector<double> delta_res;
        delta_res.push_back(p_cart[0] + dx_s_attr);
        delta_res.push_back(p_cart[1] + dy_s_attr);
        delta_res.push_back(p_cart[2] + dz_s_attr);
        Point p_res = GEOM::from_cart_cord(delta_res);
        return p_res;
    }
}

std::vector<double> ForceDirectedPP::calc_attr_spring_force_displacement(const Point &p1, const Point &p2)
{
    std::vector<double> p1_cart = GEOM::to_cart_cord(p1);
    std::vector<double> p2_cart = GEOM::to_cart_cord(p2);

    double distance = GEOM::haversine_distance_m(p1, p2);
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
    if (isnan(dx_x) || isnan(dy_y) || isnan(dz_z))
    {
        std::vector<double> p_res_cart;
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        return p_res_cart;
    }
    std::vector<double> p_res_cart;
    p_res_cart.push_back(dx_x);
    p_res_cart.push_back(dy_y);
    p_res_cart.push_back(dz_z);
    return p_res_cart;
}

std::vector<double> ForceDirectedPP::calc_rep_spring_force_displacement(const Point &p1, const Point &p2)
{
    std::vector<double> p1_cart = GEOM::to_cart_cord(p1);
    std::vector<double> p2_cart = GEOM::to_cart_cord(p2);

    double distance = GEOM::haversine_distance_m(p1, p2);
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
    std::vector<double> p1_cart = GEOM::to_cart_cord(p1);
    std::vector<double> p2_cart = GEOM::to_cart_cord(p2);
    double distance = GEOM::haversine_distance_m(p1, p2);
    // Duplicated point
    if (distance <= 0.000)
    {
        std::vector<double> p_res_cart;
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        p_res_cart.push_back(0);
        return p_res_cart;
    }
    double dx = abs(p1_cart[0] - p2_cart[0]);
    double dy = abs(p1_cart[1] - p2_cart[1]);
    double dz = abs(p1_cart[2] - p2_cart[2]);
    double fx = total_force * (dx / distance);
    double fy = total_force * (dy / distance);
    double fz = total_force * (dz / distance);
    if (total_force < 0)
    {
        if (p2_cart[0] > p1_cart[0])
            fx = -fx;
        if (p2_cart[1] > p1_cart[1])
            fy = -fy;
        if (p2_cart[2] > p1_cart[2])
            fz = -fz;
    }
    else
    {
        if (p2_cart[0] < p1_cart[0])
            fx = -fx;
        if (p2_cart[1] < p1_cart[1])
            fy = -fy;
        if (p2_cart[2] < p1_cart[2])
            fz = -fz;
    }

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
