// Generated by rstantools.  Do not edit by hand.

#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.26.1-4-gd72b68b7-dirty
#include <stan/model/model_header.hpp>
namespace model_bdSampler_namespace {
inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math;
using stan::math::pow; 
stan::math::profile_map profiles__;
static int current_statement__= 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in 'bdSampler', line 25, column 2 to column 42)",
                                                      " (in 'bdSampler', line 26, column 2 to column 33)",
                                                      " (in 'bdSampler', line 27, column 2 to column 35)",
                                                      " (in 'bdSampler', line 30, column 2 to column 44)",
                                                      " (in 'bdSampler', line 20, column 2 to column 21)",
                                                      " (in 'bdSampler', line 21, column 9 to column 14)",
                                                      " (in 'bdSampler', line 21, column 2 to column 16)",
                                                      " (in 'bdSampler', line 22, column 2 to column 19)",
                                                      " (in 'bdSampler', line 3, column 4 to column 16)",
                                                      " (in 'bdSampler', line 4, column 4 to column 12)",
                                                      " (in 'bdSampler', line 5, column 4 to column 20)",
                                                      " (in 'bdSampler', line 7, column 4 to column 11)",
                                                      " (in 'bdSampler', line 8, column 4 to column 31)",
                                                      " (in 'bdSampler', line 9, column 4 to column 88)",
                                                      " (in 'bdSampler', line 10, column 4 to column 91)",
                                                      " (in 'bdSampler', line 12, column 6 to column 23)",
                                                      " (in 'bdSampler', line 13, column 6 to column 52)",
                                                      " (in 'bdSampler', line 14, column 6 to column 95)",
                                                      " (in 'bdSampler', line 11, column 23 to line 15, column 5)",
                                                      " (in 'bdSampler', line 11, column 4 to line 15, column 5)",
                                                      " (in 'bdSampler', line 16, column 4 to column 15)",
                                                      " (in 'bdSampler', line 2, column 74 to line 17, column 3)"};
template <bool propto__, typename T0__, typename T1__, typename T2__,
typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
logLikeBDcoalTimes_lpdf(const std::vector<T0__>& t, const T1__& lambda,
                        const T2__& mu, const T3__& lgRho,
                        std::ostream* pstream__) {
  using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
  local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
  (void) DUMMY_VAR__;  // suppress unused var warning
  
  try {
    int numCoal;
    numCoal = std::numeric_limits<int>::min();
    
    local_scalar_t__ ll;
    ll = DUMMY_VAR__;
    
    current_statement__ = 11;
    numCoal = stan::math::size(t);
    current_statement__ = 12;
    ll = 0.0;
    current_statement__ = 13;
    ll = (ll + stan::math::log((numCoal + 1)));
    current_statement__ = 14;
    ll = (((ll + stan::math::log((lambda - mu))) +
            (numCoal * (stan::math::log(lambda) + lgRho))) -
           ((lambda - mu) * max(t)));
    current_statement__ = 15;
    ll = (ll -
           stan::math::log(
             ((stan::math::exp(lgRho) * lambda) +
               (((lambda * (1 - stan::math::exp(lgRho))) - mu) *
                 stan::math::exp((-(lambda - mu) * max(t)))))));
    current_statement__ = 20;
    for (int k = 1; k <= numCoal; ++k) {
      current_statement__ = 16;
      ll = (ll + stan::math::log(k));
      current_statement__ = 17;
      ll = ((ll + (2 * stan::math::log((lambda - mu)))) -
             ((lambda - mu) * t[(k - 1)]));
      current_statement__ = 18;
      ll = (ll -
             (2 *
               stan::math::log(
                 ((stan::math::exp(lgRho) * lambda) +
                   (((lambda * (1 - stan::math::exp(lgRho))) - mu) *
                     stan::math::exp((-(lambda - mu) * t[(k - 1)])))))));}
    current_statement__ = 21;
    return ll;
  } catch (const std::exception& e) {
    stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
  }
  
}
struct logLikeBDcoalTimes_lpdf_functor__ {
template <bool propto__, typename T0__, typename T1__, typename T2__,
typename T3__>
stan::promote_args_t<T0__, T1__, T2__,
T3__>
operator()(const std::vector<T0__>& t, const T1__& lambda, const T2__& mu,
           const T3__& lgRho, std::ostream* pstream__)  const 
{
return logLikeBDcoalTimes_lpdf<propto__>(t, lambda, mu, lgRho, pstream__);
}
};
#include <stan_meta_header.hpp>
class model_bdSampler final : public model_base_crtp<model_bdSampler> {
private:
  int nCoal;
  std::vector<double> t;
  double upperLambda;
 
public:
  ~model_bdSampler() { }
  
  inline std::string model_name() const final { return "model_bdSampler"; }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.26.1-4-gd72b68b7-dirty", "stancflags = "};
  }
  
  
  model_bdSampler(stan::io::var_context& context__,
                  unsigned int random_seed__ = 0,
                  std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "model_bdSampler_namespace::model_bdSampler";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      current_statement__ = 5;
      context__.validate_dims("data initialization","nCoal","int",
          context__.to_vec());
      nCoal = std::numeric_limits<int>::min();
      
      current_statement__ = 5;
      nCoal = context__.vals_i("nCoal")[(1 - 1)];
      current_statement__ = 5;
      current_statement__ = 5;
      check_greater_or_equal(function__, "nCoal", nCoal, 1);
      current_statement__ = 6;
      validate_non_negative_index("t", "nCoal", nCoal);
      current_statement__ = 7;
      context__.validate_dims("data initialization","t","double",
          context__.to_vec(nCoal));
      t = std::vector<double>(nCoal, std::numeric_limits<double>::quiet_NaN());
      
      current_statement__ = 7;
      assign(t, nil_index_list(), context__.vals_r("t"),
        "assigning variable t");
      current_statement__ = 8;
      context__.validate_dims("data initialization","upperLambda","double",
          context__.to_vec());
      upperLambda = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 8;
      upperLambda = context__.vals_r("upperLambda")[(1 - 1)];
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI, stan::require_vector_like_t<VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "model_bdSampler_namespace::log_prob";
(void) function__;  // suppress unused var warning
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      local_scalar_t__ lambda;
      lambda = DUMMY_VAR__;
      
      current_statement__ = 1;
      lambda = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        lambda = stan::math::lub_constrain(lambda, 0, upperLambda, lp__);
      } else {
        current_statement__ = 1;
        lambda = stan::math::lub_constrain(lambda, 0, upperLambda);
      }
      local_scalar_t__ mu;
      mu = DUMMY_VAR__;
      
      current_statement__ = 2;
      mu = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        mu = stan::math::lub_constrain(mu, 0, lambda, lp__);
      } else {
        current_statement__ = 2;
        mu = stan::math::lub_constrain(mu, 0, lambda);
      }
      local_scalar_t__ lgRho;
      lgRho = DUMMY_VAR__;
      
      current_statement__ = 3;
      lgRho = in__.scalar();
      current_statement__ = 3;
      if (jacobian__) {
        current_statement__ = 3;
        lgRho = stan::math::lub_constrain(lgRho, -1000, 0, lp__);
      } else {
        current_statement__ = 3;
        lgRho = stan::math::lub_constrain(lgRho, -1000, 0);
      }
      {
        current_statement__ = 4;
        lp_accum__.add(
          logLikeBDcoalTimes_lpdf<propto__>(t, lambda, mu, lgRho, pstream__));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, stan::require_std_vector_vt<std::is_floating_point, VecVar>* = nullptr>
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "model_bdSampler_namespace::write_array";
(void) function__;  // suppress unused var warning
    (void) function__;  // suppress unused var warning
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    
    try {
      double lambda;
      lambda = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      lambda = in__.scalar();
      current_statement__ = 1;
      lambda = stan::math::lub_constrain(lambda, 0, upperLambda);
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu = in__.scalar();
      current_statement__ = 2;
      mu = stan::math::lub_constrain(mu, 0, lambda);
      double lgRho;
      lgRho = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      lgRho = in__.scalar();
      current_statement__ = 3;
      lgRho = stan::math::lub_constrain(lgRho, -1000, 0);
      vars__.emplace_back(lambda);
      vars__.emplace_back(mu);
      vars__.emplace_back(lgRho);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, stan::require_std_vector_t<VecVar>* = nullptr, stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void transform_inits_impl(const stan::io::var_context& context__,
                                   VecI& params_i__, VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    vars__.clear();
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      pos__ = std::numeric_limits<int>::min();
      
      pos__ = 1;
      double lambda;
      lambda = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      lambda = context__.vals_r("lambda")[(1 - 1)];
      double lambda_free__;
      lambda_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 1;
      lambda_free__ = stan::math::lub_free(lambda, 0, upperLambda);
      double mu;
      mu = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu = context__.vals_r("mu")[(1 - 1)];
      double mu_free__;
      mu_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 2;
      mu_free__ = stan::math::lub_free(mu, 0, lambda);
      double lgRho;
      lgRho = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      lgRho = context__.vals_r("lgRho")[(1 - 1)];
      double lgRho_free__;
      lgRho_free__ = std::numeric_limits<double>::quiet_NaN();
      
      current_statement__ = 3;
      lgRho_free__ = stan::math::lub_free(lgRho, -1000, 0);
      vars__.emplace_back(lambda_free__);
      vars__.emplace_back(mu_free__);
      vars__.emplace_back(lgRho_free__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__.clear();
    names__.emplace_back("lambda");
    names__.emplace_back("mu");
    names__.emplace_back("lgRho");
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.clear();
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    dimss__.emplace_back(std::vector<size_t>{});
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "lambda");
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "lgRho");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "lambda");
    param_names__.emplace_back(std::string() + "mu");
    param_names__.emplace_back(std::string() + "lgRho");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lgRho\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lgRho\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      std::vector<double> vars_vec(vars.size());
      std::vector<int> params_i;
      write_array_impl(base_rng, params_r, params_i, vars_vec,
          emit_transformed_parameters, emit_generated_quantities, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i) {
        vars.coeffRef(i) = vars_vec[i];
      }
    }
    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      write_array_impl(base_rng, params_r, params_i, vars, emit_transformed_parameters, emit_generated_quantities, pstream);
    }
    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }
  
    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits_impl(context, params_i, params_r_vec, pstream);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i) {
        params_r.coeffRef(i) = params_r_vec[i];
      }
    }
    inline void transform_inits(const stan::io::var_context& context,
                                std::vector<int>& params_i,
                                std::vector<double>& vars,
                                std::ostream* pstream = nullptr) const final {
      transform_inits_impl(context, params_i, vars, pstream);
    }        
};
}
using stan_model = model_bdSampler_namespace::model_bdSampler;
#ifndef USING_R
// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_bdSampler_namespace::profiles__;
}
#endif
#endif
