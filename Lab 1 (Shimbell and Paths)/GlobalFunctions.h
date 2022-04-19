#pragma once
#include <string>

bool IsOnlyDigits(const std::string s)
{
	return s.find_first_not_of("0123456789") == std::string::npos;
}