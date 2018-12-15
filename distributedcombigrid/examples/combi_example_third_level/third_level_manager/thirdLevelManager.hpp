#include <iostream>
#include <string>
#include <vector>
#include <boost/property_tree/ini_parser.hpp>
#include "SimpleAmqpClient/SimpleAmqpClient.h"

int main(int argc, char* argv[]);

void readConfig();

void testMessage(std::string& queue);

void establishMessageQueues();
