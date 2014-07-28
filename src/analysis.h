#ifndef _ANALYSIS_H_
#define _ANALYSIS_H_

#include <ea/analysis.h>
#include <ea/analysis/dominant.h>
#include "sunspot.h"

LIBEA_ANALYSIS_TOOL(sunspot_data) {
    using namespace ealib;
    using namespace ealib::analysis;
    typename EA::fitness_function_type::matrix_type& M=ea.fitness_function()._train_input;
    typename EA::fitness_function_type::vector_type& T=ea.fitness_function()._train_t;
    typename EA::fitness_function_type::vector_type& Tp1=ea.fitness_function()._train_tplus1;
    
    datafile df("sunspot_data.dat");
    df.comment("WARNING: bits are reversed (i0=LSB) and multiple numbers may be encoded per line");
    df.add_field("t");
    for(std::size_t i=0; i<M.size2(); ++i) {
        df.add_field("i" + boost::lexical_cast<std::string>(i));
    }
    df.add_field("tplus1");
    
    for(std::size_t i=0; i<M.size1(); ++i) {
        typename EA::fitness_function_type::row_type r(M,i);
        df.write(T(i))
        .write_all(r.begin(), r.end())
        .write(Tp1(i))
        .endl();
    }
}


LIBEA_ANALYSIS_TOOL(sunspot_train_predictions) {
    using namespace ealib;
    using namespace ealib::analysis;
    namespace bnu=boost::numeric::ublas;
    
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_train_predictions.dat");
    
    df.add_field("sample").add_field("t").add_field("observed").add_field("predicted");
    
    sunspot_fitness::matrix_type output;
    ea.fitness_function().train(ind, output, ea.rng(), ea);
    
    sunspot_fitness::vector_type& observed = ea.fitness_function()._test_tplus1;
    double factor = static_cast<double>(1 << get<SUNSPOT_FRACTIONAL_BITS>(ea));
    
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        sunspot_fitness::column_type c(output,i);
        assert(c.size() == observed.size());
        
        sunspot_fitness::dvector_type obs = bnu::vector_range<sunspot_fitness::vector_type>(observed, bnu::range(i,observed.size())) / factor;
        sunspot_fitness::dvector_type pre = bnu::vector_range<sunspot_fitness::column_type>(c, bnu::range(0,c.size()-i)) / factor;
        assert(obs.size() == pre.size());
        
        for(std::size_t j=0; j<obs.size(); ++j) {
            df.write(j).write(i+1).write(obs(j)).write(pre(j)).endl();
        }
    }
}

LIBEA_ANALYSIS_TOOL(sunspot_test_predictions) {
    using namespace ealib;
    using namespace ealib::analysis;
    namespace bnu=boost::numeric::ublas;
    
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_test_predictions.dat");
    
    df.add_field("sample").add_field("t").add_field("observed").add_field("predicted");
    
    sunspot_fitness::matrix_type output;
    ea.fitness_function().test(ind, output, ea.rng(), ea);
    
    sunspot_fitness::vector_type& observed = ea.fitness_function()._test_tplus1;
    double factor = static_cast<double>(1 << get<SUNSPOT_FRACTIONAL_BITS>(ea));
    
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        sunspot_fitness::column_type c(output,i);
        assert(c.size() == observed.size());
        
        sunspot_fitness::dvector_type obs = bnu::vector_range<sunspot_fitness::vector_type>(observed, bnu::range(i,observed.size())) / factor;
        sunspot_fitness::dvector_type pre = bnu::vector_range<sunspot_fitness::column_type>(c, bnu::range(0,c.size()-i)) / factor;
        assert(obs.size() == pre.size());
        
        for(std::size_t j=0; j<obs.size(); ++j) {
            df.write(j).write(i+1).write(obs(j)).write(pre(j)).endl();
        }
    }
}


LIBEA_ANALYSIS_TOOL(sunspot_train_detail) {
    using namespace ealib;
    using namespace ealib::analysis;
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_train_detail.dat");
    df.add_field("observed");
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
    }
    
    sunspot_fitness::matrix_type output;
    ea.fitness_function().train(ind, output, ea.rng(), ea);
    
    for(std::size_t i=0; i<output.size1(); ++i) {
        df.write(ea.fitness_function()._train_tplus1(i));
        
        for(std::size_t j=0; j<output.size2(); ++j) {
            df.write(output(i,j));
        }
        
        df.endl();
    }
}


LIBEA_ANALYSIS_TOOL(sunspot_test_detail) {
    using namespace ealib;
    using namespace ealib::analysis;
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_test_detail.dat");
    df.add_field("observed");
    
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
    }
    
    sunspot_fitness::matrix_type output;
    ea.fitness_function().test(ind, output, ea.rng(), ea);
    
    for(std::size_t i=0; i<output.size1(); ++i) {
        df.write(ea.fitness_function()._test_tplus1(i));
        
        for(std::size_t j=0; j<output.size2(); ++j) {
            df.write(output(i,j));
        }
        
        df.endl();
    }
}

LIBEA_ANALYSIS_TOOL(sunspot_train_rmse) {
    using namespace ealib;
    using namespace ealib::analysis;
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_train_rmse.dat");
    df.add_field("total_rmse");
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
    }
    
    sunspot_fitness::matrix_type output;
    sunspot_fitness::dvector_type rmse = ea.fitness_function().train(ind, output, ea.rng(), ea);
    
    df.write(std::accumulate(rmse.begin(),rmse.end(),0.0)).write_all(rmse.begin(), rmse.end()).endl();
}


LIBEA_ANALYSIS_TOOL(sunspot_test_rmse) {
    using namespace ealib;
    using namespace ealib::analysis;
    typename EA::individual_type& ind = *dominant(ea);
    
    datafile df("sunspot_test_rmse.dat");
    df.add_field("total_rmse");
    for(std::size_t i=0; i<get<SUNSPOT_PREDICTION_HORIZON>(ea); ++i) {
        df.add_field("tplus" + boost::lexical_cast<std::string>(i+1));
    }
    
    sunspot_fitness::matrix_type output;
    sunspot_fitness::dvector_type rmse = ea.fitness_function().test(ind, output, ea.rng(), ea);
    
    df.write(std::accumulate(rmse.begin(),rmse.end(),0.0)).write_all(rmse.begin(), rmse.end()).endl();
}

#endif
