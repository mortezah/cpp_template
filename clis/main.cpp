#include <iostream>
#include <format>
#include <ranges>

#include "vector_algebra_utils.h"

int main(int argc, char** argv)
{
  std::cout << std::format("number is {:04} and text is {}\n", 31, "hellos!");
  return 0;
}
