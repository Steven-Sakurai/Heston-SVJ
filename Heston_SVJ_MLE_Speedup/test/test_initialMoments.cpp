#include "InitialMoments.hpp"

int main() {
	stdVec M(4*21);
	stdVec my_par = {0.05, 2, 0.2, 0.25, -0.8};
	initialMoments(my_par, M, true);
	return 0;
}
