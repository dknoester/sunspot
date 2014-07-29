/* sunspot.cpp
 *
 * This file is part of sunspot.
 *
 * Copyright 2014 David B. Knoester.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <ea/mkv/markov_network_evolution.h>
#include <ea/generational_models/wright_fisher.h>
#include <ea/fitness_function.h>
#include <ea/cmdline_interface.h>
#include <ea/datafiles/fitness.h>
using namespace ealib;
using namespace mkv;

#include "sunspot.h"
#include "analysis.h"


typedef markov_network_evolution
< sunspot_fitness
, recombination::asexual
, generational_models::wright_fisher< >
> ea_type;


template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void gather_options() {
        add_mkv_options(this);
        
        // ea options
        add_option<POPULATION_SIZE>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        // sunspot options
        add_option<SUNSPOT_TRAIN>(this);
        add_option<SUNSPOT_TEST>(this);
        add_option<SUNSPOT_INTEGER_BITS>(this);
        add_option<SUNSPOT_FRACTIONAL_BITS>(this);
        add_option<SUNSPOT_PREDICTION_HORIZON>(this);
        add_option<SUNSPOT_LIMIT>(this);
        add_option<SUNSPOT_INPUT_LAGS>(this);
        add_option<SUNSPOT_ENCODING>(this);
    }
    
    virtual void gather_tools() {
        add_tool<sunspot_data>(this);
        add_tool<sunspot_train_detail>(this);
        add_tool<sunspot_test_detail>(this);
        add_tool<sunspot_train_rmse>(this);
        add_tool<sunspot_test_rmse>(this);
        add_tool<sunspot_train_predictions>(this);
        add_tool<sunspot_test_predictions>(this);
        add_tool<analysis::dominant_genetic_graph>(this);
        add_tool<analysis::dominant_causal_graph>(this);
        add_tool<analysis::dominant_reduced_graph>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<datafiles::fitness_dat>(ea);
    }

    virtual void after_initialization(EA& ea) {
        std::size_t nbits = get<SUNSPOT_INTEGER_BITS>(ea) + get<SUNSPOT_FRACTIONAL_BITS>(ea);
        assert(nbits < (sizeof(long)*8));
        assert(get<SUNSPOT_INPUT_LAGS>(ea) >= 1);

        put<MKV_INPUT_N>(get<SUNSPOT_INPUT_LAGS>(ea) * nbits, ea);
        put<MKV_OUTPUT_N>(2 * nbits * get<SUNSPOT_PREDICTION_HORIZON>(ea), ea);
    }
};
LIBEA_CMDLINE_INSTANCE(ea_type, cli);
