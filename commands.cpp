
#include "global.h"
#include "supercool.h"
#include "cpptcl/cpptcl.h"

#include "tcc7A.h"
#include "tcc9A.h"
#include "tcc10B.h"
#include "tcc11A.h"
#include "tcc11W.h"
#include "tcc13A.h"
#include "tccBCC9.h"

#include "trajectory_func.h"
#include "wrap_particles.h"

/**
 * Determine clusters from static TCC. Select particles bound in clusters.
 */
class tf_clusters_tcc : public trajectory_func
{
private:
  std::string name;
  bool is;
  bool select_special;

protected:
  virtual double apply_cb( dynamics *dyn, particle_vector &pv )
  {
    tcc_base *tcc;
    if( name == "7A" )
      tcc = new tcc7A( pv.size() );
    else if( name == "9A" )
      tcc = new tcc9A( pv.size() );
    else if( name == "10B" )
      tcc = new tcc10B( pv.size() );
    else if( name == "11A" )
      tcc = new tcc11A( pv.size() );
    else if( name == "11W" )
      tcc = new tcc11W( pv.size() );
    else if( name == "13A" )
      tcc = new tcc13A( pv.size() );
    else if( name == "BCC9" )
      tcc = new tccBCC9( pv.size() );
    else {
      // error: unknown cluster type
      std::string msg = "error: unknown cluster type ";
      msg += name;
      g_ctx.error( msg.c_str() );
      return 0.0;
    }

    tcc->set_particles( dyn->get_box(), pv, is );
    tcc->Bonds_GetBondsV();
    int N = tcc->Find();
    // set cluster info
    pv.clusters = tcc->clusters;
    // legacy code: set empty info for every cluster
    if( pv.clusters.empty() )
      pv.clusters = std::vector<clusterinfo>( tcc->nClu );
    // mark particles in clusters
    for( int n = 0; n < pv.size(); n++ ) {
      pv[n].data = 0.0;
      if( select_special )
        pv[n].selected = tcc->sSpe[n] ? true : false;
      else
        pv[n].selected = tcc->sClu[n] ? true : false;
    }
    delete tcc;

    return N;
  }

public:
  tf_clusters_tcc( const char *_name, const char *_rng, int flags )
    : name( _name )
  {
    rng.compile( _rng );
    is = flags & 1;
    select_special = flags & 2;
  }

  virtual double operator()( dynamics *dyn, trajectory *traj )
  {
    return apply( dyn, traj );
  }

  static trajectory_func *create( const char *name, const char *rng, int flags )
  {
    return new tf_clusters_tcc( name, rng, flags );
  }
};

// ==== compatibility layer ====

trajectory_func *tcc9A( const char *rng, double rc, bool is )
{
  return new tf_clusters_tcc( "9A", rng, is );
}

trajectory_func *tcc11A( const char *rng, double rc, bool is )
{
  return new tf_clusters_tcc( "11A", rng, is );
}

trajectory_func *tcc13A( const char *rng, double rc, bool is )
{
  return new tf_clusters_tcc( "13A", rng, is );
}

trajectory_func *tccBCC9( const char *rng, double rc, bool is )
{
  return new tf_clusters_tcc( "BCC9", rng, is );
}

// ==== expose ====

void expose_commands_clusters( Tcl::interpreter &i )
{
  using namespace Tcl;

  i.def( "tf_clusters_tcc", &tf_clusters_tcc::create, factory( "trajectory_func" ) );

  i.def( "tf_clusters_tcc9A", &tcc9A, factory( "trajectory_func" ) );
  i.def( "tf_clusters_tcc11A", &tcc11A, factory( "trajectory_func" ) );
  i.def( "tf_clusters_tcc13A", &tcc13A, factory( "trajectory_func" ) );
  i.def( "tf_clusters_tccBCC9", &tccBCC9, factory( "trajectory_func" ) );
}
