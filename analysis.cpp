#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;

// Среднеквадратичное отклонение...

double calculateMean(const std::vector<double>& data) {
    double sum = 0.0;
    for (double value : data) {
        sum += value;
    }
    return sum / data.size();
}

double calculateStandardDeviation(const std::vector<double>& data) {
    double mean = calculateMean(data);
    double sumSquaredDifferences = 0.0;
    for (double value : data) {
        sumSquaredDifferences += std::pow(value - mean, 2);
    }
    double variance = sumSquaredDifferences / data.size();
    return std::sqrt(variance);
}

void simu( int argc, char **argv )
{
    std::ifstream is_geom( "geom1.dat"  );
    if( !is_geom.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[1] + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();
    geom.build_surface();

    std::ifstream is_epot( "epot1.dat" );
    if( !is_epot.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[2] + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    std::ifstream is_pdb( "pdb1.dat" );
    if( !is_pdb.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[3] + "\'" ) );
    ParticleDataBase3D pdb( is_pdb, geom );
    is_pdb.close();

    VectorField *bfield = NULL;
	// bool fout[3] = {true, true, true};
	// MeshVectorField *mesh_bfield = new MeshVectorField( MODE_3D, fout, 1.0e-3, 1.0, "bfield1.txt");
	// field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO, 
	// 				 FIELD_ZERO, FIELD_ZERO, 
	// 				 FIELD_ZERO, FIELD_ZERO };
	// mesh_bfield->set_extrapolation( bfldextrpl );
	// bfield = mesh_bfield;
    int n_step = 500;
    double step = ( geom.max(2) - geom.origo(2) )/n_step ;
    TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_X );
        diagnostics.push_back( DIAG_XP );
        diagnostics.push_back( DIAG_Y );
        diagnostics.push_back( DIAG_YP );

    
    // for( size_t k = 0; k < n_step+1; k++ ) {
    //     pdb.trajectories_at_plane( tdata, AXIS_Z, geom.origo(2) + step*k, diagnostics );
    //     Emittance emitx( tdata(0).data(), tdata(1).data() );
    //     Emittance emity( tdata(2).data(), tdata(3).data() );
    //     // Output
        
    //     int vec_size = tdata(0).data().size();   

    //     std::cout<<tdata(0).data().size(); 

    //     double rms_x = calculateStandardDeviation(tdata(0).data());
    //     double mean_x = calculateMean(tdata(0).data());
    //     double rms_y = calculateStandardDeviation(tdata(2).data());
    //     double mean_y = calculateMean(tdata(2).data());

    //     // std::cout << "RMS для положительных значений: " << rms_positive*2 << std::endl;
    //     // std::cout << "RMS для отрицательных значений: " << -rms_negative*2 << std::endl;

    //     // ofstream dout1( "x_z_noscc.txt", ios_base::app );
    //     //     dout1 << geom.origo(2) + step*k << " "
    //     //     << mean_x + rms_x*2 << " "
    //     //     << mean_x - rms_x*2 << "\n";

    //     // ofstream dout2( "y_z_noscc.txt", ios_base::app );
    //     //     dout2 << geom.origo(2) + step*k << " "
    //     //     << mean_y + rms_y*2 << " "
    //     //     << mean_y - rms_y*2 << "\n";

    //     // ofstream dout3( "emit_noscc.txt", ios_base::app );
    //     //     dout3 << geom.origo(2) + step*k << " "
    //     //     << emitx.epsilon() << " "
    //     //     << emity.epsilon() << "\n";

    //     // ofstream dout( "y_z.txt", ios_base::app );
    //     //     dout << geom.origo(2) + step*k << " "
    //     //     << tdata(0).data()[i] << "\n";

    //     // for(int i=0;i<vec_size+1;i++){
    //     //     ofstream dout( "x_z.txt", ios_base::app );
    //     //     dout << geom.origo(2) + step*k << " "
    //     //     << tdata(0).data()[i] << "\n";
    //     // }

    // }


       

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
    if( bfield )
	plotter.set_bfield( bfield );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    // if( argc <= 3 ) {
	// cerr << "Usage: analysis geom epot pdb <bfield>\n";
	// exit( 1 );
    // }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
