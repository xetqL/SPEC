//
// Created by xetql on 5/6/19.
//

#ifndef SPEC_CLIPARSER_HPP
#define SPEC_CLIPARSER_HPP

#include <zupply.hpp>
#include <SimulationParams.hpp>
#include <utility>
std::pair<SimulationParams, bool>  parse_cli(int argc, char** argv);

#endif //SPEC_CLIPARSER_HPP
