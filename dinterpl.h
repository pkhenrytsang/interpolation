/***

    Interpolation Routine by Pak Ki Henry Tsang
    Use at your own risk!

***/

#include <cstdlib>
#include <cstring>

using namespace std;

class dinterpl{

  public:
  
    /* evaluation cache */
    typedef struct {
    size_t  cache;        /* cache of index   */
    size_t  cache_hit;   /* keep statistics  */
    size_t  cache_miss;
    }cache;
    
    static cache * cache_alloc (void)
    {
      cache *a = (cache *) malloc (sizeof (cache));

      a->cache = 0;
      a->cache_hit = 0;
      a->cache_miss = 0;

      return a;
    }
    static void cache_free (cache * a)
    {
      free (a);
    }

    //Class constructor/destructor
    dinterpl(const double x_array[], const double y_array[], size_t size);
    ~dinterpl();
    
    //Interpolation Evaluation
    double linear_eval(double x, cache * cache);
    double cspline_eval(double x, cache * cache);
    
    double linear_eval(double x);
    double cspline_eval(double x);
    
    //Misc Routines
    void uniformize_grid(size_t s_size);
    

  private:
  
    //Arrays stored in memory (must free them with destructor)
    double* cspline_a;
    double* cspline_b;
    double* cspline_c;
    double* cspline_d;
    double* x_array;
    double* y_array;

    size_t size;

    size_t internal_cache;
    
    void cspline_init(const double x_array[], const double y_array[], size_t size);


    /*
      Inline search fuctions
    */
    inline static size_t interp_bsearch (const double x_array[], double x, size_t index_lo, size_t index_hi) 
    {
      size_t min = index_lo, max = index_hi;
      while (min + 1 < max)
      {
        size_t i = (min + max) >> 1;
        min = x > x_array[i] ? i : min;
        max = x > x_array[i] ? max : i;
      }
      return min;
    }

    inline static size_t interp_accel_find(cache * cache, const double xa[], size_t len, double x)
    {
      size_t x_index = cache->cache;

      if(x < xa[x_index]) {
        cache->cache = interp_bsearch(xa, x, 0, x_index);
        cache->cache_miss++;
      }
      else if(x >= xa[x_index + 1]) {
        cache->cache = interp_bsearch(xa, x, x_index, len-1);
        cache->cache_miss++;
      }
      else {
        cache->cache_hit++;
      }

      return cache->cache;
    }

    inline static size_t interp_accel_find(size_t& internal_cache, const double xa[], size_t len, double x)
    {
      size_t x_index = internal_cache;

      if(x < xa[x_index]) {
        internal_cache = interp_bsearch(xa, x, 0, x_index);
      }
      else if(x >= xa[x_index + 1]) {
        internal_cache = interp_bsearch(xa, x, x_index, len-1);
      }

      return internal_cache;
    }
    
    inline static double cspline_eval_at_knot(size_t i, double t
                                              ,const double cspline_a[],const double cspline_b[]
                                              ,const double cspline_c[],const double cspline_d[])
    {
      return cspline_a[i]+cspline_b[i]*t+cspline_c[i]*t*t+cspline_d[i]*t*t*t;
    }
    
    void free_memory();

};

