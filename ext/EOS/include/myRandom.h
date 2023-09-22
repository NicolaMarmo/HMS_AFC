#ifndef RANDOM_MERSENNE
#define RANDOM_MERSENNE

void mysgenrand(unsigned long seed);
double mygenrand();

//#include<random>
//int intRand(const int & min, const int & max) {
//	static thread_local std::mt19937 generator;
//	std::uniform_int_distribution<int> distribution(min, max);
//	return distribution(generator);
//}

#endif // RANDOM_MERSENNE