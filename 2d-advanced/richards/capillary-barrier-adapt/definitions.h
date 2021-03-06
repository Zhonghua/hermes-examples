#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Constitutive relations.
enum CONSTITUTIVE_RELATIONS {
    CONSTITUTIVE_GENUCHTEN,    // Van Genuchten.
    CONSTITUTIVE_GARDNER       // Gardner.
};

// Problem parameters.
extern double H_INIT;                            // Initial pressure head.
extern double H_ELEVATION;                        // Top constant pressure head -- an infiltration experiment.
extern double K_S_vals[4]; 
extern double ALPHA_vals[4];
extern double N_vals[4];
extern double M_vals[4];

extern double THETA_R_vals[4];
extern double THETA_S_vals[4];
extern double STORATIVITY_vals[4];


// Precalculation of constitutive tables.
extern const int MATERIAL_COUNT;
extern const int CONSTITUTIVE_TABLE_METHOD;          // 0 - constitutive functions are evaluated directly (slow performance).
					                                        // 1 - constitutive functions are linearly approximated on interval <TABLE_LIMIT; LOW_LIMIT> 
						                                      //	 (very efficient CPU utilization less efficient memory consumption (depending on TABLE_PRECISION)).
						                                      // 2 - constitutive functions are aproximated by quintic splines.
						  
//!Use only if 	CONSTITUTIVE_TABLE_METHOD == 2 !//					  
extern const int NUM_OF_INTERVALS;                  // Number of intervals.                      
extern const double INTERVALS_4_APPROX[16]; 
//!------------------------------------------!//


//!Use only if CONSTITUTIVE_TABLE_METHOD == 1 !//
extern double TABLE_LIMIT; 		                // Limit of precalculated functions (should be always negative value lower 
						                                      // then the lowest expect value of the solution (consider DMP!!)
extern const double TABLE_PRECISION;               // Precision of precalculated table use 1.0, 0,1, 0.01, etc.....
extern const double LOW_LIMIT;                      // Lower bound of K(h) function approximated by polynomials.
extern const int NUM_OF_INSIDE_PTS;
//!------------------------------------------!//


extern bool POLYNOMIALS_READY;
extern bool POLYNOMIALS_ALLOCATED;

class ConstitutiveRelations
{
  public:
  ConstitutiveRelations(double table_limit);

  // This function uses Horner scheme to efficiently evaluate values of 
  // polynomials in approximating K(h) function close to full saturation.
  virtual double K(double h, int layer) = 0;
  
  virtual double dKdh(double h, int layer) = 0;
  
  virtual double ddKdhh(double h, int layer) = 0;
  
  virtual double C(double h, int layer) = 0;
  
  virtual double dCdh(double h, int layer) = 0;

  virtual bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer) = 0;

  // Members.
  bool polynomials_ready;
  bool constitutive_tables_ready;
  double table_limit;
};

class ConstitutiveGenuchten : public ConstitutiveRelations
{
public:
  ConstitutiveGenuchten(double low_limit, bool polynomials_ready, int constitutive_table_method, int num_of_inside_pts, 
                        double table_limit, double* alpha_vals, double* n_vals, double* m_vals, double* k_s_vals, 
                        double* theta_r_vals, double* theta_s_vals, double* storativity_vals, double table_precision, int material_count, bool polynomials_allocated,
                        int num_of_intervals, int iterative_method): ConstitutiveRelations(table_limit),
   low_limit(low_limit), polynomials_ready(polynomials_ready), constitutive_table_method(constitutive_table_method), polynomials(polynomials), 
     num_of_inside_pts(num_of_inside_pts), alpha_vals(alpha_vals), 
     n_vals(n_vals), m_vals(m_vals), k_s_vals(k_s_vals), theta_r_vals(theta_r_vals), theta_s_vals(theta_s_vals), storativity_vals(storativity_vals),  
     table_precision(table_precision), material_count(material_count), polynomials_allocated(polynomials_allocated), num_of_intervals(num_of_intervals)
  {
    constitutive_tables_ready = false;
    if(constitutive_table_method == 1) {
      this->get_constitutive_tables(iterative_method);
      constitutive_tables_ready = true;
    }
  }

  // This function uses Horner scheme to efficiently evaluate values of 
  // polynomials in approximating K(h) function close to full saturation.
  double horner(double *pol, double x, int n){
    double px=0.0;
    for (int i=0; i<n; i++){
      px = px*x + pol[n-1-i];
    }
    return px;
  }

  // K (van Genuchten).
  virtual double K(double h, int layer) {
    double value ;
    int location ;
    if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1){
      return horner(polynomials[layer][0], h, (6+num_of_inside_pts));
    }
    else {
      if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit)  {
        alpha = alpha_vals[layer];
        n = n_vals[layer];
        m = m_vals[layer];
        k_s = k_s_vals[layer];
        theta_r = theta_r_vals[layer];
        theta_s = theta_s_vals[layer];
        storativity = storativity_vals[layer];
        if (h < 0) return 
	    k_s*(std::pow(1 - std::pow(-(alpha*h), m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m),2)/
        std::pow(1 + std::pow(-(alpha*h), n), m/2.));
        else return k_s;    
      }
      else if (h<0 && constitutive_table_method == 1) {
        location = -int(h/table_precision) ;
        value = (K_table[layer][location+1] - K_table[layer][location])*(-h/table_precision-location)+K_table[layer][location];
        return value;
      }
      else if (h<0 && constitutive_table_method ==2) {
        location = pol_search_help[int(-h)]; 
        return horner(k_pols[location][layer][0], h, 6);
      }
      else return k_s_vals[layer];
    }
  }

  // dK/dh (van Genuchten).
  virtual double dKdh(double h, int layer) {
    double value ;
    int location ;
  
    if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1){
      return horner(polynomials[layer][1], h, (5+num_of_inside_pts));
    }
    else {
      if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit)  {
        alpha = alpha_vals[layer];
        n = n_vals[layer];
        m = m_vals[layer];
        k_s = k_s_vals[layer];
        theta_r = theta_r_vals[layer];
        theta_s = theta_s_vals[layer];
        storativity = storativity_vals[layer];
        if (h < 0) return 
	    k_s*((alpha*std::pow(-(alpha*h), -1 + n)*
	    std::pow(1 + std::pow(-(alpha*h), n), -1 - m/2.)*
	    std::pow(1 - std::pow(-(alpha*h), m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m), 2)*m*n)/2. + 
	    (2*(1 - std::pow(-(alpha*h), m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m))*
	    (-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
	    std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) + 
	    (alpha*std::pow(-(alpha*h), -1 + m*n)*m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m)))/
	    std::pow(1 + std::pow(-(alpha*h), n), m/2.));
        else return 0;
      }
      else if (h < 0 && constitutive_table_method == 1) {
        location = -int(h/table_precision);
        value = (dKdh_table[layer][location+1] - dKdh_table[layer][location])*(-h/table_precision-location)+dKdh_table[layer][location] ;
        return value ;
      }
      else if (h<0 && constitutive_table_method ==2) {
        location = pol_search_help[int(-h)]; 
        return horner(k_pols[location][layer][1], h, 5);
      }
      else return 0 ;
    }
  }

  // ddK/dhh (van Genuchten).
  virtual double ddKdhh(double h, int layer) {

    int location ;
    double value ;
  
    if (h > low_limit && h < 0 && polynomials_ready && constitutive_table_method == 1){
      return horner(polynomials[layer][2], h, (4+num_of_inside_pts));
    }
    else {
      if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
        alpha = alpha_vals[layer];
        n = n_vals[layer];
        m = m_vals[layer];
        k_s = k_s_vals[layer];
        theta_r = theta_r_vals[layer];
        theta_s = theta_s_vals[layer];
        storativity = storativity_vals[layer];

        if (h  < 0 ) return 
	      k_s*( -(std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n)*
	    std::pow(1 + std::pow(-(alpha*h), n), -1 - m/2.)*
	    std::pow(1 - std::pow(-(alpha*h), m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m), 2)*m*(-1 + n)*n)/2.\
	      - (std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2*n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -2 - m/2.)*
	    std::pow(1 - std::pow(-(alpha*h), m*n)/
	        std::pow(1 + std::pow(-(alpha*h), n), m), 2)*(-1 - m/2.)*m*
	    std::pow(n, 2))/2. + 2*alpha*std::pow(-(alpha*h), -1 + n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -1 - m/2.)*
	    (1 - std::pow(-(alpha*h), m*n)/
	      std::pow(1 + std::pow(-(alpha*h), n), m))*m*n*
	    (-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) + 
	    (alpha*std::pow(-(alpha*h), -1 + m*n)*m*n)/
	    std::pow(1 + std::pow(-(alpha*h), n), m)) + 
	    (2*std::pow(-(alpha*std::pow(-(alpha*h), -1 + n + m*n)*
	        std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n) + 
	      (alpha*std::pow(-(alpha*h), -1 + m*n)*m*n)/
	      std::pow(1 + std::pow(-(alpha*h), n), m), 2))/
	      std::pow(1 + std::pow(-(alpha*h), n), m/2.) + 
	    (2*(1 - std::pow(-(alpha*h), m*n)/
	      std::pow(1 + std::pow(-(alpha*h), n), m))*
	      (std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2*n + m*n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -2 - m)*(-1 - m)*m*
	      std::pow(n, 2) + std::pow(alpha, 2)*
	      std::pow(-(alpha*h), -2 + n + m*n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*std::pow(m, 2)*
	      std::pow(n, 2) - (std::pow(alpha, 2)*
	        std::pow(-(alpha*h), -2 + m*n)*m*n*(-1 + m*n))/
	      std::pow(1 + std::pow(-(alpha*h), n), m) + 
	      std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n + m*n)*
	      std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*
	      (-1 + n + m*n)))/std::pow(1 + std::pow(-(alpha*h), n), m/2.));

        else return 0;
      }
      else if (h < 0 && constitutive_table_method == 1) {
        location = -int(h/table_precision) ;
        value = (ddKdhh_table[layer][location+1] - 
                ddKdhh_table[layer][location])*(-h/table_precision-location)+ddKdhh_table[layer][location];
        return value ;
      }
      else if (h<0 && constitutive_table_method == 2) {
        location = pol_search_help[int(-h)]; 
        return horner(k_pols[location][layer][2], h, 4);
      }
      else return 0;
    }	
  }

  // C (van Genuchten).
  virtual double C(double h, int layer)
  {
    int location ;
    double value ;
  
    if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit)  {
      alpha = alpha_vals[layer];
      n = n_vals[layer];
      m = m_vals[layer];
      k_s = k_s_vals[layer];
      theta_r = theta_r_vals[layer];
      theta_s = theta_s_vals[layer];
      storativity = storativity_vals[layer];
      if (h < 0) return 
      alpha*std::pow(-(alpha*h), -1 + n)*
        std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*(theta_s - theta_r) + 
      (storativity*((theta_s - theta_r)/std::pow(1 + std::pow(-(alpha*h), n), m) + 
	    theta_r))/theta_s;
      else return storativity; 
    }
    else if (h<0 && constitutive_table_method == 1) {
      location = -int(h/table_precision);
      value = (C_table[layer][location+1] - C_table[layer][location])*(-h/table_precision-location)+C_table[layer][location];
      return value;
    }
    else if (h<0 && constitutive_table_method == 2) {
        location = pol_search_help[int(-h)]; 
        return horner(c_pols[location][layer][0], h, 4);
    }
    else return storativity_vals[layer];
  }

  // dC/dh (van Genuchten).
  virtual double dCdh(double h, int layer)
  {
    int location ;
    double value ;
    
    if (constitutive_table_method == 0 || !constitutive_tables_ready || h < table_limit) {
      alpha = alpha_vals[layer];
      n = n_vals[layer];
      m = m_vals[layer];
      k_s = k_s_vals[layer];
      theta_r = theta_r_vals[layer];
      theta_s = theta_s_vals[layer];
      storativity = storativity_vals[layer];
      if (h < 0) return
              -(std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + n)*
	        std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*(-1 + n)*n*
	        (theta_s - theta_r)) - std::pow(alpha, 2)*std::pow(-(alpha*h), -2 + 2*n)*
              std::pow(1 + std::pow(-(alpha*h), n), -2 - m)*(-1 - m)*m*std::pow(n, 2)*
              (theta_s - theta_r) + (alpha*std::pow(-(alpha*h), -1 + n)*
	        std::pow(1 + std::pow(-(alpha*h), n), -1 - m)*m*n*storativity*
	        (theta_s - theta_r))/theta_s;
      else return 0; 
    }
    else if (h<0 && constitutive_table_method == 1) {
      location = -int(h/table_precision) ;
      value = (dCdh_table[layer][location+1] - 
              dCdh_table[layer][location])*(-h/table_precision-location)+dCdh_table[layer][location];
      return value ;
    }
    else if (h<0 && constitutive_table_method == 2) {
        location = pol_search_help[int(-h)]; 
        return horner(c_pols[location][layer][1], h, 3);
    }
    else return 0;
      
  }

  // Initialize polynomial approximation of constitutive relations close to full saturation for CONSTITUTIVE_TABLE_METHOD=1.
  // For CONSTITUTIVE_TABLE_METHOD=2 all constitutive functions are approximated by polynomials, K(h) function by quintic spline, C(h) 
  // function by cubic splines. Discretization is managed by variable int NUM_OF_INTERVALS and double* INTERVALS_4_APPROX.
  // ------------------------------------------------------------------------------
  // For CONSTITUTIVE_TABLE_METHOD=1 this function requires folowing arguments:
  // n - degree of polynomials
  // low_limit - start point of the polynomial approximation
  // points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
  // An approriate amount of points related to the polynomial degree should be selected.
  // n_inside_point - number of inside points
  // layer - material to be considered.
  //------------------------------------------------------------------------------
  // For CONSTITUTIVE_TABLE_METHOD=2, all parameters are obtained from global definitions.
  virtual bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer) {
    double** Aside;
    double* Bside;
    switch (constitutive_table_method) { 
      // no approximation 
      case 0 :
	      break ;
      // polynomial approximation only for the the K(h) function surroundings close to zero
      case 1 :
	      if (polynomials_allocated == false) {
	        polynomials = new double**[material_count];
	        for (int i=0; i<material_count; i++)
	          polynomials[i] = new double*[3];
	        for (int i=0; i<material_count; i++)
	          for (int j=0; j<3; j++)
	            polynomials[i][j] = new double[n_inside_points+6] ;
	        polynomials_allocated = true;
	      }

	      Aside = new double*[n+n_inside_points];
	      Bside = new double[n+n_inside_points];
	      for (int i=0; i<n; i++)
	        Aside[i] = new double[n+n_inside_points];

	      // Evaluate the first three rows of the matrix (zero, first and second derivative at point low_limit).
	      for (int i=0; i<n; i++)
	        for (int j=0; j<n; j++)
	          Aside[i][j] = 0.0;
	      for (int i=0; i<n; i++){
	        Aside[3][i] = std::pow(low_limit,i) ;
	        Aside[4][i] = i*std::pow(low_limit, i-1) ;
	        Aside[5][i] = i*(i-1)*std::pow(low_limit, i-2) ;
	      }
	      Bside[3] = K(low_limit, layer) ;
	      Bside[4] = dKdh(low_limit, layer) ;
	      Bside[5] = ddKdhh(low_limit, layer) ; 
	
	      // Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	      Aside[0][0] = 1.0 ;

	      // For the both first and second derivative it does not really matter what value is placed there.
	      Aside[1][1] = 1.0 ;
	      Aside[2][2] = 2.0 ;
      
	      Bside[0] = K(0.0, layer) ;
	      Bside[1] = 0.0;
	      Bside[2] = 0.0;
      
	      for (int i=6; i<(6+n_inside_points); i++) {
	        for (int j=0; j<n; j++)
	          Aside[i][j]=std::pow(points[i-6],j) ;
	        
          printf("poradi, %i %lf %lf \n", i, K(points[i-6], layer), points[i-6]);

	        Bside[i] = K(points[i-6], layer) ;
	        printf("layer, %i \n", layer);
	      }

	      gem_full(Aside, Bside, polynomials[layer][0], (n_inside_points+6));
	
	      for (int i=1; i<3; i++){
	        for (int j=0; j< (n_inside_points+5); j++)
	          polynomials[layer][i][j] = (j+1)*polynomials[layer][i-1][j+1];
	        polynomials[layer][i][n_inside_points+6-i] = 0.0;
	      }
	
	      delete [] Aside;
	      delete [] Bside;
	      break ;
      // polynomial approximation for all functions at interval (table_limit, 0)
      case 2 :
	      int pts = 0;
	      if (polynomials_allocated == false) {
	        // K(h) function is approximated by quintic spline.
	        k_pols = new double***[num_of_intervals];
	        //C(h) function is approximated by cubic spline.
	        c_pols = new double***[num_of_intervals];
	  
	        for (int i=0; i<num_of_intervals; i++) {
	          k_pols[i] = new double**[material_count];
	          c_pols[i] = new double**[material_count];
	          for (int j=0; j<material_count; j++) {
	            k_pols[i][j] = new double*[3];
	            c_pols[i][j] = new double*[2];
	            for (int k=0; k<3; k++) {
		            k_pols[i][j][k] = new double[6 + pts];
		            if (k<2)
		              c_pols[i][j][k] = new double[4];
	            }
	          }
	        }
	  
  	      //allocate pol_search_help array -- an index array with locations for particular pressure head functions
	        pol_search_help = new int[int(-table_limit)+1];
	  
	        for (int i=0; i<=int(-table_limit); i++) {
	          for (int j=0; j<num_of_intervals; j++) {
	            if (j < 1) {
		            if (-i > INTERVALS_4_APPROX[j]) {
		              pol_search_help[i] = j ;
		              break ;
		            }
	            }
	            else {
		            if (-i >= INTERVALS_4_APPROX[j] && -i <= INTERVALS_4_APPROX[j-1]) {
		              pol_search_help[i] = j ;
		              break ;
		            }
	            }
	          }
	        }
	  
	        polynomials_allocated = true;
	      }
	      //create matrix
	      Aside = new double*[6 + pts];
	      for (int i=0; i< (6 + pts); i++)
	        Aside[i] = new double[6+pts];
	
	      Bside = new double[6 + pts];
        for (int i=0; i<num_of_intervals; i++) {
	        if (i < 1) {
	          for (int j=0; j<3; j++)
	            for (int k=0; k<(6 + pts); k++)
		            Aside[j][k] = 0.0;
	         // Evaluate the second three rows of the matrix (zero, first and second derivative at point zero).
	          Aside[0][0] = 1.0 ;

	          // For the both first and second derivative it does not really matter what value is placed there.
	          Aside[1][1] = 1.0 ;
	          Aside[2][2] = 2.0 ;
	  
	        }
	        else {
	          for (int j=0; j<(6+pts); j++) {
	            Aside[0][j] = std::pow(INTERVALS_4_APPROX[i-1],j);
	            Aside[1][j] = j*std::pow(INTERVALS_4_APPROX[i-1], j-1);
	            Aside[2][j] = j*(j-1)*std::pow(INTERVALS_4_APPROX[i-1], j-2);
	          }
	        }
	        for (int j=0; j<(6+pts); j++) {
	          Aside[3][j] = std::pow(INTERVALS_4_APPROX[i],j);
	          Aside[4][j] = j*std::pow(INTERVALS_4_APPROX[i], j-1);
	          Aside[5][j] =  j*(j-1)*std::pow(INTERVALS_4_APPROX[i], j-2);
	          switch (pts) {
	            case 0 :
		            break;
	            case 1 : 
		            if (i > 0)
		              Aside[6][j] =  std::pow((INTERVALS_4_APPROX[i]+INTERVALS_4_APPROX[i-1])/2, j);
		            else
		              Aside[6][j] =  std::pow((INTERVALS_4_APPROX[i])/2, j) ;
		            break;
	            default :
		            printf("too many of inside points in polynomial approximation; not implemented!!! (check extras.cpp) \n");
		            exit(1);
	          }
	        }
	        //Evaluate K(h) function.
	        if (i<1) {
	          Bside[0] = K(0.0, layer);
	          Bside[1] = 0.0;
	          Bside[2] = 0.0;
	        }
	        else {
	          Bside[0] = K(INTERVALS_4_APPROX[i-1], layer);
	          Bside[1] = dKdh(INTERVALS_4_APPROX[i-1], layer);
	          Bside[2] = ddKdhh(INTERVALS_4_APPROX[i-1], layer);
	        }
	 
	        Bside[3] = K(INTERVALS_4_APPROX[i], layer);
	        Bside[4] = dKdh(INTERVALS_4_APPROX[i], layer);
	        Bside[5] = ddKdhh(INTERVALS_4_APPROX[i], layer);  
	  
	        switch (pts) {
	            case 0 :
		            break;
	            case 1 :
		            if (i > 0)
		              Bside[6] = K((INTERVALS_4_APPROX[i]+INTERVALS_4_APPROX[i-1])/2, layer);
		            else
		              Bside[6] = K((INTERVALS_4_APPROX[i])/2, layer);
		            break;
	        }

	        gem_full(Aside, Bside, k_pols[i][layer][0], (6+pts));

	        for (int j=1; j<3; j++) {
	          for (int k=0; k<5; k++)
	            k_pols[i][layer][j][k] = (k+1)*k_pols[i][layer][j-1][k+1];
	          k_pols[i][layer][j][6-j] = 0.0;
	        }
	  
	        //Evaluate C(h) functions.
	        if (i<1) {
	          Bside[0] = C(0.0, layer);
	          Bside[1] = dCdh(0.0, layer);
	        }
	        else {
	          Bside[0] = C(INTERVALS_4_APPROX[i-1], layer);
	          Bside[1] = dCdh(INTERVALS_4_APPROX[i-1], layer);
	        }
	 
	        //The first two lines of the matrix Aside stays the same.
	        for (int j=0; j<4; j++) {	   
	          Aside[2][j] = std::pow(INTERVALS_4_APPROX[i],j);
	          Aside[3][j] = j*std::pow(INTERVALS_4_APPROX[i], j-1);
	        }

	  
	        Bside[2] = C(INTERVALS_4_APPROX[i], layer);
	        Bside[3] = dCdh(INTERVALS_4_APPROX[i], layer);
	 
	  
	        gem_full(Aside, Bside, c_pols[i][layer][0], 4);
	  
	        for (int k=0; k<5; k++)
	          c_pols[i][layer][1][k] = (k+1)*c_pols[i][layer][0][k+1];
	  
	        c_pols[i][layer][1][5] = 0.0;
	
        }
      delete [] Aside;
      delete [] Bside;
      break ;
    }
    return true;
  }

  bool get_constitutive_tables(int method)
  {
    info("Creating tables of constitutive functions (complicated real exponent relations).");

    // Table values dimension.
    int bound = int(-table_limit/table_precision)+1;
  
    // Allocating arrays. 
    K_table = new double*[material_count];
    for (int i=0; i<material_count; i++)
      K_table[i] = new double[bound];
  
    dKdh_table = new double*[material_count] ;
    for (int i=0; i<material_count; i++)
      dKdh_table[i] = new double[bound];
  
    dKdh_table = new double*[material_count] ;
    for (int i=0; i<material_count; i++)
      dKdh_table[i] = new double[bound];

    C_table = new double*[material_count] ;
    for (int i=0; i<material_count; i++)
      C_table[i] = new double[bound];
  
    //If Newton method (method==1) selected constitutive function derivations are required.
    if (method==1){
      dCdh_table = new double*[material_count] ;
      for (int i=0; i<material_count; i++)
        dCdh_table[i] = new double[bound];
    
      ddKdhh_table = new double*[material_count] ;
      for (int i=0; i<material_count; i++)
        ddKdhh_table[i] = new double[bound];
    }
  
    // Calculate and save K(h).
    info("Calculating and saving K(h).");
    for (int j=0; j<material_count; j++)
      for (int i=0; i< bound; i++)
        K_table[j][i] = K(-table_precision*i, j);
    // Calculate and save dKdh(h).
    info("Calculating and saving dKdh(h).");
    for (int j=0; j<material_count; j++)
      for (int i=0; i< bound; i++)
        dKdh_table[j][i] = dKdh(-table_precision*i, j);
    
    // Calculate and save C(h).
    info("Calculating and saving C(h).");
    for (int j=0; j<material_count; j++)
      for (int i=0; i< bound; i++)
        C_table[j][i] = C(-table_precision*i, j);
  
    //If Newton method (method==1) selected constitutive function derivations are required.
    if (method==1){
      // Calculate and save ddKdhh(h).
      info("Calculating and saving ddKdhh(h).");
      for (int j=0; j<material_count; j++)
        for (int i=0; i< bound; i++)
	        ddKdhh_table[j][i] = ddKdhh(-table_precision*i, j);
      // Calculate and save dCdh(h).
      info("Calculating and saving dCdh(h).");
      for (int j=0; j<material_count; j++)
        for (int i=0; i< bound; i++)
	        dCdh_table[j][i] = dCdh(-table_precision*i, j);
    }	
      
    return true;
  }

  // Simple Gaussian elimination for full matrices called from init_polynomials().
  bool gem_full(double** A, double* b, double* X, int n) {
    int i,j,k;
  
    double** aa;
    double dotproduct, tmp;
    aa = new double*[n];
  
    for (i=0; i<n; i++){
      aa[i] = new double[n+1];
    }
   
    for (i=0;i<n; i++){
      for (j=0; j<n; j++){
        aa[i][j] = A[i][j] ;
      }
      aa[i][n] = b[i];
    }

    for (j=0; j<(n-1); j++){
      for (i=j+1; i<n; i++){
      tmp = aa[i][j]/aa[j][j];

        for (k=0; k<(n+1); k++){
	  aa[i][k] = aa[i][k] - tmp*aa[j][k] ;
        }
      }
    }

    for (i=n-1; i>-1; i--){
      dotproduct=0.0;
      for (j=i+1; j<n; j++){
        dotproduct = dotproduct + aa[i][j]*X[j] ;
      }
      X[i] = (aa[i][n]-dotproduct)/aa[i][i] ;
    }
    
    delete []aa;
    return true;
  }

  // members.
  double low_limit; 
  bool polynomials_ready; 
  int constitutive_table_method; 
  double*** polynomials;  // Polynomial approximation of the K(h) function close to saturation.
                          // This function has singularity in its second derivative.
						              // First dimension is material ID
						              // Second dimension is the polynomial derivative.
						              // Third dimension are the polynomial's coefficients.
  int num_of_inside_pts; 
  double* alpha_vals; 
  double* n_vals; 
  double* m_vals; 
  double* k_s_vals; 
  double* theta_r_vals; 
  double* theta_s_vals; 
  double* storativity_vals; 
  double table_precision;
  int material_count;
  bool polynomials_allocated;
  int num_of_intervals;
  double alpha;
  double n;
  double m;
  double k_s;
  double theta_r;
  double theta_s;
  double storativity;
  double** K_table;
  double** dKdh_table;
  double** ddKdhh_table;
  double** C_table;
  double** dCdh_table;
  double**** c_pols;
  double**** k_pols;
  int* pol_search_help;

  bool constitutive_tables_ready;
};

class ConstitutiveGardner : public ConstitutiveRelations
{
public:
  ConstitutiveGardner(double k_s, double alpha, double theta_s, double theta_r, int constitutive_table_method) : ConstitutiveRelations(0.0),
      k_s(k_s), alpha(alpha), theta_s(theta_s), theta_r(theta_r), constitutive_table_method(constitutive_table_method) {}
  
  // K (Gardner).
  virtual double K(double h, int layer) {
    if (h < 0) return k_s*exp(alpha*h);
    else return k_s;    
  }

  // dK/dh (Gardner).
  virtual double dKdh(double h, int layer) {
    if (h < 0) return k_s*alpha*exp(alpha*h);
    else return 0;
  }

  // ddK/dhh (Gardner).
  virtual double ddKdhh(double h, int layer) {
    if (h < 0) return k_s*alpha*alpha*exp(alpha*h);
    else return 0;
  }

  // C (Gardner).
  virtual double C(double h, int layer) {
    if (h < 0) return alpha*(theta_s - theta_r)*exp(alpha*h);
    else return alpha*(theta_s - theta_r); 
  //   else return storativity; 
  }

  // dC/dh (Gardner).
  virtual double dCdh(double h, int layer) {
    if (h < 0) return alpha*(theta_s - theta_r)*alpha*exp(alpha*h);
    else return 0;    
  }

  // Initialize polynomial approximation of constitutive relations close to full saturation for CONSTITUTIVE_TABLE_METHOD=1.
  // For CONSTITUTIVE_TABLE_METHOD=2 all constitutive functions are approximated by polynomials, K(h) function by quintic spline, C(h) 
  // function by cubic splines. Discretization is managed by variable int NUM_OF_INTERVALS and double* INTERVALS_4_APPROX.
  // ------------------------------------------------------------------------------
  // For CONSTITUTIVE_TABLE_METHOD=1 this function requires folowing arguments:
  // n - degree of polynomials
  // low_limit - start point of the polynomial approximation
  // points - array of points inside the interval bounded by <low_limit, 0> to improve the accuracy, at least one is recommended. 
  // An approriate amount of points related to the polynomial degree should be selected.
  // n_inside_point - number of inside points
  // layer - material to be considered.
  //------------------------------------------------------------------------------
  // For CONSTITUTIVE_TABLE_METHOD=2, all parameters are obtained from global definitions.
  virtual bool init_polynomials(int n, double low_limit, double *points, int n_inside_points, int layer) {
    switch (constitutive_table_method) { 
      // no approximation 
      case 0 :
	      break;
      // polynomial approximation only for the the K(h) function surroundings close to zero
      case 1 :
	      exit(1);
      case 2 :
	      exit(2);
    }
    return true;
  }

  // Members.
  double k_s;
  double alpha;
  double theta_s;
  double theta_r;
  int constitutive_table_method;
};

class ConstitutiveGenuchten;
class ConstitutiveGardner;

// The first part of the file dontains forms for the Newton's
// method. Forms for Picard are in the second part.

/*** INITIAL CONDITION ***/

class InitialSolutionRichards : public ExactSolutionScalar<double>
{
public:
  InitialSolutionRichards(Mesh* mesh, double constant) 
         : ExactSolutionScalar<double>(mesh), constant(constant) {};

  virtual double value (double x, double y) const {
    return constant;
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 0.0;
    dy = 0.0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(0);
  }

  // Value.
  double constant;
};

class ExactSolutionPoisson : public ExactSolutionScalar<double>
{
public:
  ExactSolutionPoisson(Mesh* mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual double value (double x, double y) const {
    return x*x +y*y;
  }

  virtual void derivatives (double x, double y, double& dx, double& dy) const {
    dx = 2*x;
    dy = 2*y;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return x*x +y*y;
  }
};

/*** NEWTON ***/

class WeakFormRichardsNewtonEuler : public WeakForm<double>
{
public:
  WeakFormRichardsNewtonEuler(ConstitutiveRelations* relations, double tau, Solution<double>* prev_time_sln, Mesh* mesh) 
               : WeakForm<double>(1), mesh(mesh) {
    JacobianFormNewtonEuler* jac_form = new JacobianFormNewtonEuler(0, 0, relations, tau);
    jac_form->ext.push_back(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonEuler* res_form = new ResidualFormNewtonEuler(0, relations, tau);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonEuler : public MatrixFormVol<double>
  {
  public:
    JacobianFormNewtonEuler(int i, int j, ConstitutiveRelations* relations, double tau) 
      : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
                         std::string elem_marker = static_cast<WeakFormRichardsNewtonEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];

      for (int i = 0; i < n; i++)
        result += wt[i] * (
        relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                               * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                          );
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormNewtonEuler : public VectorFormVol<double>
  {
  public:
    ResidualFormNewtonEuler(int i, ConstitutiveRelations* relations, double tau) 
               : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsNewtonEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * (
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) 
                           * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           -relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          );
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  Mesh* mesh;
};


class WeakFormRichardsNewtonCrankNicolson : public WeakForm<double>
{
public:
  WeakFormRichardsNewtonCrankNicolson(ConstitutiveRelations* relations, double tau, Solution<double>* prev_time_sln, Mesh* mesh) : WeakForm<double>(1), mesh(mesh) {
    JacobianFormNewtonCrankNicolson* jac_form = new JacobianFormNewtonCrankNicolson(0, 0, relations, tau);
    jac_form->ext.push_back(prev_time_sln);
    add_matrix_form(jac_form);

    ResidualFormNewtonCrankNicolson* res_form = new ResidualFormNewtonCrankNicolson(0, relations, tau);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormNewtonCrankNicolson : public MatrixFormVol<double>
  {
  public:
    JacobianFormNewtonCrankNicolson(int i, int j, ConstitutiveRelations* relations, double tau) 
      : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsNewtonCrankNicolson*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++)
        result += wt[i] * 0.5 * ( // implicit Euler part:
		             relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
		             + relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->val[i] * v->val[i] / tau
		             - relations->dCdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_time->val[i] * v->val[i] / tau
			     + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             + relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * 
                               (h_prev_newton->dx[i]*v->dx[i] + h_prev_newton->dy[i]*v->dy[i])
                             - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]
                             - relations->ddKdhh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * u->val[i] * h_prev_newton->dy[i] * v->val[i]
                           )
                + wt[i] * 0.5 * ( // explicit Euler part, 
		             relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                           );
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormNewtonCrankNicolson : public VectorFormVol<double>
  {
  public:
    ResidualFormNewtonCrankNicolson(int i, ConstitutiveRelations* relations, double tau) : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      double result = 0;
      std::string elem_marker = static_cast<WeakFormRichardsNewtonCrankNicolson*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      Func<double>* h_prev_newton = u_ext[0];
      Func<double>* h_prev_time = ext->fn[0];
      for (int i = 0; i < n; i++) {
        result += wt[i] * 0.5 * ( // implicit Euler part
		           relations->C(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_newton->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->dx[i] * v->dx[i] + h_prev_newton->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_newton->val[i], atoi(elem_marker.c_str())) * h_prev_newton->dy[i] * v->val[i]
                          )
                + wt[i] * 0.5 * ( // explicit Euler part
		           relations->C(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_newton->val[i] - h_prev_time->val[i]) * v->val[i] / tau
                           + relations->K(h_prev_time->val[i], atoi(elem_marker.c_str())) * (h_prev_time->dx[i] * v->dx[i] + h_prev_time->dy[i] * v->dy[i])
                           - relations->dKdh(h_prev_time->val[i], atoi(elem_marker.c_str())) * h_prev_time->dy[i] * v->val[i]
		           );
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  Mesh* mesh;
};


class WeakFormRichardsPicardEuler : public WeakForm<double>
{
public:
  WeakFormRichardsPicardEuler(ConstitutiveRelations* relations, double tau, Solution<double>* prev_picard_sln, Solution<double>* prev_time_sln, Mesh* mesh) : WeakForm<double>(1), mesh(mesh) {
    JacobianFormPicardEuler* jac_form = new JacobianFormPicardEuler(0, 0, relations, tau);
    jac_form->ext.push_back(prev_picard_sln);
    add_matrix_form(jac_form);

    ResidualFormPicardEuler* res_form = new ResidualFormPicardEuler(0, relations, tau);
    res_form->ext.push_back(prev_picard_sln);
    res_form->ext.push_back(prev_time_sln);
    add_vector_form(res_form);
  }

private:
  class JacobianFormPicardEuler : public MatrixFormVol<double>
  {
  public:
    JacobianFormPicardEuler(int i, int j, ConstitutiveRelations* relations, double tau) 
      : MatrixFormVol<double>(i, j, HERMES_ANY, HERMES_NONSYM), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsPicardEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_picard = ext->fn[0];

      for (int i = 0; i < n; i++) {
        result += wt[i] * (  relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->val[i] * v->val[i] / tau
                             + relations->K(h_prev_picard->val[i], atoi(elem_marker.c_str())) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                             - relations->dKdh(h_prev_picard->val[i], atoi(elem_marker.c_str())) * u->dy[i] * v->val[i]);
      }
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return matrix_form<double, double>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }

    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  class ResidualFormPicardEuler : public VectorFormVol<double>
  {
  public:
    ResidualFormPicardEuler(int i, ConstitutiveRelations* relations, double tau) : VectorFormVol<double>(i), tau(tau), relations(relations) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      std::string elem_marker = static_cast<WeakFormRichardsPicardEuler*>(wf)->mesh->get_element_markers_conversion().get_user_marker(e->elem_marker).marker;
      double result = 0;
      Func<double>* h_prev_picard = ext->fn[0];
      Func<double>* h_prev_time = ext->fn[1];
      for (int i = 0; i < n; i++) 
        result += wt[i] * relations->C(h_prev_picard->val[i], atoi(elem_marker.c_str())) * h_prev_time->val[i] * v->val[i] / tau;
      return result;
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<double> *ext) const {
      return vector_form<double, double>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return Ord(30);
    }
    
    // Members.
    double tau;
    ConstitutiveRelations* relations;
  };

  Mesh* mesh;
};

class RichardsEssentialBC : public EssentialBoundaryCondition<double> {
public:

  RichardsEssentialBC(std::string marker, double h_elevation, double pulse_end_time, double h_init, double startup_time) :
  EssentialBoundaryCondition<double>(Hermes::vector<std::string>()), h_elevation(h_elevation), pulse_end_time(pulse_end_time), h_init(h_init), startup_time(startup_time)
  {
    markers.push_back(marker);
  }

  ~RichardsEssentialBC() {}

  inline EssentialBCValueType get_value_type() const { return EssentialBoundaryCondition<double>::BC_FUNCTION; }

  virtual double value(double x, double y, double n_x, double n_y, double t_x, double t_y) const {
    if (current_time < startup_time)
      return h_init + current_time/startup_time*(h_elevation-h_init);
    else if (current_time > pulse_end_time)
      return h_init;
    else
      return h_elevation;
  }

  // Member.
  double h_elevation;
  double pulse_end_time;
  double h_init;
  double startup_time;
};