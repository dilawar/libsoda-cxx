#ifndef HELPER_H
#define HELPER_H

/* --------------------------------------------------------------------------*/
/**
 * @Synopsis  Check the equality of two numbers.
 *
 * @Param a
 * @Param b
 *
 * @Returns   
 */
/* ----------------------------------------------------------------------------*/
template<typename T=double>
bool areEqual(T a, T b) 
{
    return std::fabs(a - b) < 1e-6;
}



#endif /* end of include guard: HELPER_H */
