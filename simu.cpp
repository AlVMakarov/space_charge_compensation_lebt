#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "multimeshvectorfield.hpp"
#include "meshvectorfield.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include "readascii.hpp"
#include "scharge.hpp"
#include <ctime>

using namespace std;

/*
Инициализация магнитных полей
*/
const std::string bfieldfn_sol1 = "SOL1.txt";
const std::string bfieldfn_sol2 = "SOL2.txt";

/*
Геометрия канала - апертура x мм
*/
bool solid(double x, double y, double z)
{
    double r = sqrt(x * x + y * y);
    return (r >= 0.065);
}
/*
Обновление статуса частицы (для корректного расчета динамики электронов)
*/
void reset_p_stat(ParticleDataBase3D &pdb)
{
    int k;
    for (int i = 0; i < pdb.size(); i++)
    {
        k = pdb.traj_size(i);
        Particle3D &p = pdb.particle(i);
        if (p.get_status() != PARTICLE_OK)
        {
            p.set_status(PARTICLE_OK);
        }
    };
}

/*
Потенциал z=0, phi(x)
*/
void plot_epot(const std::vector<double> &x, const std::vector<double> &y, const std::string &picname)
{

    // Открываем файл для записи данных
    std::ofstream data_file("data.dat");

    // Записываем данные в файл
    for (size_t i = 0; i < x.size(); ++i)
    {
        data_file << x[i] << " " << y[i] << std::endl;
    }

    // Закрываем файл
    data_file.close();

    // Формируем команду для gnuplot
    string command = "gnuplot -e \"set terminal png; set output '" + picname + "'; plot 'data.dat' with lines\"";

    // Выполняем команду gnuplot
    system(command.c_str());

    // Удаляем файл с данными
    std::remove("data.dat");
}

/*
Построение графика плотности траекторий
*/
void plot_trajectory_density(const Geometry &geom, const ParticleDataBase3D &pdb, EpotField &epot,
                             const MeshScalarField scharge, const std::string &picname)
{
    MeshScalarField tdens(geom);
    pdb.build_trajectory_density_field(tdens);
    GeomPlotter gplotter(geom);
    gplotter.set_size(4096, 2048);
    gplotter.set_view(VIEW_ZX);
    gplotter.set_epot(&epot);
    gplotter.set_font_size(10);
    gplotter.set_scharge(&scharge);
    std::vector<double> eqlines;
    eqlines.push_back(-4.0);
    eqlines.push_back(-2.0);
    eqlines.push_back(0.01);
    eqlines.push_back(+2.0);
    eqlines.push_back(+4.0);
    gplotter.set_eqlines_manual(eqlines);
    gplotter.set_particle_database(&pdb);
    gplotter.set_particle_div(0);
    gplotter.set_trajdens(&tdens);
    // gplotter.set_fieldgraph_plot(FIELD_TRAJDENS);
    gplotter.set_fieldgraph_plot(FIELD_SCHARGE);
    gplotter.fieldgraph()->set_zscale(ZSCALE_RELLOG);
    gplotter.plot_png(picname);
}

/*
Построение графика траекторий
*/
void plot_trajectory(const Geometry &geom, const ParticleDataBase3D &pdb, EpotField &epot, const std::string &picname)
{
    MeshScalarField tdens(geom);
    pdb.build_trajectory_density_field(tdens);
    GeomPlotter gplotter(geom);
    gplotter.set_size(4096, 2048);
    gplotter.set_view(VIEW_ZX);
    gplotter.set_epot(&epot);
    gplotter.set_font_size(10);
    std::vector<double> eqlines;
    eqlines.push_back(-4.0);
    eqlines.push_back(-2.0);
    eqlines.push_back(0.01);
    eqlines.push_back(+2.0);
    eqlines.push_back(+4.0);
    gplotter.set_eqlines_manual(eqlines);
    gplotter.set_particle_database(&pdb);
    gplotter.set_particle_div(1);
    gplotter.plot_png(picname);
}

void simu(int *argc, char ***argv)
{
    srand((int)time(0));

    double h = 5e-4;
    Vec3D origo(-70e-3,
                -70e-3,
                0.e-3);
    double sizereq[3] = {140e-3,
                         140e-3,
                         300e-3};
    Int3D meshsize((int)floor(sizereq[0] / h) + 1,
                   (int)floor(sizereq[1] / h) + 1,
                   (int)floor(sizereq[2] / h) + 1);

    Geometry geom(MODE_3D, meshsize, origo, h);

    Solid *s1 = new FuncSolid(solid);

    geom.set_solid(7, s1);

    geom.set_boundary(1, Bound(BOUND_NEUMANN, 0.0));
    geom.set_boundary(2, Bound(BOUND_NEUMANN, 0.0));
    geom.set_boundary(3, Bound(BOUND_NEUMANN, 0.0));
    geom.set_boundary(4, Bound(BOUND_NEUMANN, 0.0));
    geom.set_boundary(5, Bound(BOUND_NEUMANN, 0.0));
    geom.set_boundary(6, Bound(BOUND_NEUMANN, 0.0));

    geom.set_boundary(7, Bound(BOUND_DIRICHLET, 0.0));

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver(geom);
    EpotField epot(geom);                   // Сетка потенциалов
    MeshScalarField scharge(geom);          // Сетка зарядов ионов
    MeshScalarField scharge_ave(geom);      // Усредн сетка
    MeshScalarField scharge_ele(geom);      // Сетка зарядов электронов
    MeshScalarField scharge_ele_buff(geom); // Сетка суммарного заряда

    /*
    Инициализациия магнитных полей
    */
    bool fout[3] = {true, true, true};
    field_extrpl_e bfldextrpl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                    FIELD_ZERO, FIELD_ZERO};

    // double bfield_step = 0.001;
    // Int3D meshsize_bfield((int)floor(sizereq[0] / bfield_step) + 1,
    //                       (int)floor(sizereq[1] / bfield_step) + 1,
    //                       (int)floor(sizereq[2] / bfield_step) + 1);

    // MeshVectorField bfield(MODE_3D, fout, meshsize_bfield, origo, bfield_step);
    // MeshVectorField bfield1n(MODE_3D, fout, 1.0e-3, 16.0 / 16.0, bfieldfn_sol1);
    // MeshVectorField bfield2n(MODE_3D, fout, 1.0e-3, 16.0 / 16.0, bfieldfn_sol2);

    // bfield.set_extrapolation(bfldextrpl);
    // bfield1n.set_extrapolation(bfldextrpl);
    // bfield2n.set_extrapolation(bfldextrpl);

    // MeshVectorField bfield1(MODE_3D, fout, meshsize_bfield, origo, bfield_step, bfield1n);
    // MeshVectorField bfield2(MODE_3D, fout, meshsize_bfield, origo, bfield_step, bfield2n);

    // bfield += bfield1;
    // bfield += bfield2;
    // bfield1.clear();
    // bfield2.clear();
    MeshVectorField bfield;

    /*
    Инициализация сетки электрических полей
    */
    EpotEfield efield(epot);
    EpotEfield efield_test(epot);
    field_extrpl_e efldextrpl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                    FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE};
    efield.set_extrapolation(efldextrpl);
    efield_test.set_extrapolation(efldextrpl);

    /*
    Инициализация дб частиц
    */

    // Ионы
    ParticleDataBase3D pdb(geom);
    bool pmirror[6] = {false, false, false, false, false, false};
    pdb.set_mirror(pmirror);
    pdb.set_surface_collision(true);

    // Электроны
    ParticleDataBase3D pdb_elec(geom);
    pdb_elec.set_mirror(pmirror);
    pdb_elec.set_surface_collision(true);

    /*
    Параметры пучка, на входе в канал LEBT
    */
    const int N_COUNT = 200; // Число частиц
    // const double M = 1.0;              // Масса частицы H
    const double M = 1.0;
    const double Q = 1.0;    // Заряд частицы
    const double E0 = 25000; // Энергия частицы, эВ
    double R0 = 20e-3;       // Радиус пучка, м
    // 2.5 A - ток пучка, при котором достигается потенциал пучка в центре -> 30 В
    double I = 14e-3;                  // Ток пучка, А
    double J_R = I / (R0 * R0 * M_PI); // Плотность тока, А/м^2

    /*
    Кусок кода, отвечающий за добавление вторичных электронов в симуляцию
    */
    double meanFreePath = 1;
    double elecMass = 1.0 / 1836.00;
    double dt = 1.0e-7;
    double simuTime = 2.0e-7;
    double currTime = 0.0;
    double iTime = 0.0;

    double a;
    double elecL;
    double elecZ;
    char picname[20];

    for (int i = 0; i < 2; i++)
    {
        solver.solve(epot, scharge);
        efield.recalculate();

        int iternum = 0;
        char picname1[20];

        /*
        Начальное приближение поля
        */
        pdb.clear();
        pdb.add_cylindrical_beam_with_energy(N_COUNT, J_R, Q, M,
                                             E0, 0, 0,
                                             Vec3D(0, 0, geom.origo(2)),
                                             Vec3D(1, 0, 0),
                                             Vec3D(0, 1, 0), R0);

        pdb.iterate_trajectories(scharge, efield, bfield);

        // Построение графиков
        sprintf(picname1, "trajdens_i_check%d.png", i);
        plot_trajectory_density(geom, pdb, epot, scharge, picname1);

        double z0t = 0.1;
        double y0t = 0;
        double x0t = origo[0];
        double x1t = x0t + sizereq[0];
        std::vector<double> xArrayt;
        std::vector<double> epArrayt;
        while (x0t <= x1t + h)
        {
            xArrayt.push_back(x0t);
            epArrayt.push_back(epot(Vec3D(x0t, y0t, z0t)));
            x0t += h;
        }
        // Построение графиков
        sprintf(picname, "phi_test_%d.png", i);
        plot_epot(xArrayt, epArrayt, picname);
        
        // Уср сетки зарядов
        // if (i == 1)
        // {
        //     scharge_ave = scharge;
        // }
        // else
        // {
        //     uint32_t nodecnt = scharge.nodecount();
        //     for (uint32_t b = 0; b < nodecnt; b++)
        //     {
        //         scharge_ave(b) = 0.5 * scharge(b) + 0.5 * scharge_ave(b);
        //     }
        // }
    }

    // solver.solve(epot, scharge);
    // efield.recalculate();

    // int iter;
    /*
    Динамика вторичных частиц
    Делаем n количество итераций с заданным шагом
    */
    //int iter = 0;
    //while (false)
     for (int iter = 0; iter < 101; iter++)
    {
        // Логи
        ibsimu.message(1) << "General cycle, iter:  " << iter << "\n";
        ibsimu.message(1) << "-----------------------\n";
        solver.solve(epot, scharge);
        efield.recalculate();
        std::cout << "\n---ELec count: " << pdb_elec.size() << "--- \n";
        

        pdb.clear();
        pdb.add_cylindrical_beam_with_energy(N_COUNT, J_R, Q, M,
                                             E0, 0, 0,
                                             Vec3D(0, 0, geom.origo(2)),
                                             Vec3D(1, 0, 0),
                                             Vec3D(0, 1, 0), R0);

        // pdb.iterate_trajectories(scharge, efield, bfield);

        if (iter == 0)
        {
            pdb.iterate_trajectories(scharge, efield, bfield);
        }
        else
        {
            pdb.iterate_trajectories(scharge, efield, bfield);
        }

        for (int i = 0; i < N_COUNT; i++)
        {
            a = (double)(rand()) / RAND_MAX;
            elecL = meanFreePath * a;
            int trajSize = pdb.traj_size(i); // Coun of i-th traj point
            elecZ = elecL;                   // elec path z
            Vec3D loc0(0, 0, 0), loc1(0, 0, 0);
            Vec3D vel0(0, 0, 0), vel1(0, 0, 0);
            double time0, time1;

            for (int j = 1; j < trajSize; j++)
            {

                pdb.trajectory_point(time0, loc0, vel0, i, j);
                pdb.trajectory_point(time1, loc1, vel1, i, j-1);
                double iQCurr = pdb.particle(i).IQ()*1.0e-11;
                //double iQCurr = pdb.particle(i).IQ();
                std::cout<<iQCurr<<"\n" ;               
                //double iQCurr = 1.6e-19;
                double xCurr = loc0[0];
                double yCurr = loc0[1];
                double zCurr = loc0[2];
                double timeCurr = iter * dt;
                while (elecZ < zCurr)
                {
                    double el_energy = 0; // eV
                    double cosphi = (-1. + 2.*(double)(rand()) / RAND_MAX);
                    double sinphi = sqrt(1 - pow(cosphi,2));
                    double costet = (-1. + 2.*(double)(rand()) / RAND_MAX);
                    double xV = sqrt(2.0 * el_energy * CHARGE_E / (elecMass * MASS_U)) * cosphi;
                    double yV = sqrt(2.0 * el_energy * CHARGE_E / (elecMass * MASS_U)) * sinphi;
                    double zV = sqrt(2.0 * el_energy * CHARGE_E / (elecMass * MASS_U)) * costet;

                    // double iQCurr = scharge(Vec3D(xCurr, yCurr, zCurr))/1000;
                    // std::cout<<iQCurr<<" fgwegwe\n";
                    pdb_elec.add_particle(iQCurr, -1.0, elecMass, ParticleP3D(0, xCurr, xV, yCurr, yV, elecZ, zV));
                    // pdb_elec.add_particle(iQCurr, -1.0, elecMass, ParticleP3D(0.0, xCurr, 0.0, yCurr, 0.0, elecZ, 0.0));
                    //pdb_elec.add_particle(iQCurr, -1.0, elecMass, ParticleP3D(0.0, xCurr, 0.0, yCurr, 0.0, elecZ, 0.0));

                    a = (double)(rand()) / RAND_MAX;
                    elecL = meanFreePath * a;
                    elecZ = elecZ + elecL;
                }
            }
        }

        solver.solve(epot, scharge);
        efield.recalculate();
        /*
        Динамика электронов
        */
        // pdb_elec.clear_trajectories();
        // reset_p_stat(pdb_elec);
        // pdb_elec.set_max_time(dt * (iter + 1));
        // pdb_elec.iterate_trajectories(scharge_ele, efield, bfield);
        scharge_ele.clear();
        //pdb_elec.clear_trajectories();
        while (currTime < dt * (iter + 1))
        {
            pdb_elec.step_particles(scharge_ele_buff, efield, bfield, 1.0e-11);
            scharge_ele += scharge_ele_buff;
            currTime += 1.0e-11;
            //std::cout<<"Done \n";
        }
        scharge_finalize_pic(scharge);

        // /*
        // Дебаг
        // */
        // for (int j = 0; j < pdb_elec.size(); j++)
        // {
        //     Particle3D &pp = pdb_elec.particle(j);
        //     pp.debug_print(std::cout);
        // }

        // Построение графиков
        sprintf(picname, "trajdens_e%d.png", iter);
        plot_trajectory_density(geom, pdb_elec, epot, scharge_ele, picname);
        sprintf(picname, "traj_e%d.png", iter);
        plot_trajectory(geom, pdb_elec, epot, picname);
        sprintf(picname, "trajdens_i%d.png", iter);
        plot_trajectory_density(geom, pdb, epot, scharge, picname);
        sprintf(picname, "traj_i%d.png", iter);
        plot_trajectory(geom, pdb, epot, picname);
        // Суммирование SC
        // uint32_t nodecount = scharge.nodecount();
        // for (uint32_t b = 0; b < nodecount; b++)
        // {
        //     scharge(b) = scharge(b) + scharge_ele(b);
        // }
        scharge += scharge_ele;
        plot_trajectory_density(geom, pdb, epot, scharge, "sum_pot.png");
        // Потенциал в z=0.3, phi(x)
        // Графики
        double z0 = 0.6;
        double y0 = 0;
        double x0 = origo[0];
        double x1 = x0 + sizereq[0];
        std::vector<double> xArray;
        std::vector<double> epArray;
        while (x0 <= x1 + h)
        {
            xArray.push_back(x0);
            epArray.push_back(epot(Vec3D(x0, y0, z0)));
            x0 += h;
        }
        // Построение графиков
        sprintf(picname, "phi_%d.png", iter);
        plot_epot(xArray, epArray, picname);
    }

    MeshScalarField tdens(geom);
    pdb.build_trajectory_density_field(tdens);
    GTKPlotter plotter(argc, argv);
    plotter.set_geometry(&geom);
    plotter.set_epot(&epot);
    plotter.set_bfield(&bfield);
    plotter.set_trajdens(&tdens);
    plotter.set_scharge(&scharge);
    // plotter.set_scharge( &scharge_ele   );
    plotter.set_particledatabase(&pdb);
    // plotter.set_particledatabase(&pdb_elec);
    plotter.new_geometry_plot_window();
    plotter.run();
}

int main(int argc, char **argv)
{

    try
    {
        ibsimu.set_message_threshold(MSG_VERBOSE, 1);
        ibsimu.set_thread_count(10);
        simu(&argc, &argv);
    }
    catch (Error e)
    {
        e.print_error_message(ibsimu.message(0));
        exit(1);
    }

    return (0);
}
