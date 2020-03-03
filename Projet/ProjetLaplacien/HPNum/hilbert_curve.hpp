#ifndef _HILBERT_CURVE_HPP_
#define _HILBERT_CURVE_HPP_
#include <stdexcept>
#include <iostream>
#include <vector>

class HilbertCurve
{
public:
    HilbertCurve( unsigned niter, unsigned dim ) : number_of_iterations(niter), dimension(dim)
    {
        this->max_h = (1ULL<<(niter*dim))-1;
        this->max_x = (1ULL<<niter)-1;
    }

    /**
     * @brief      Return the coordinates for a given hilbert distance.
     *
     * @param[in]  h     integer distance along hilbert curve
     *
     * @return     transpose of h (n components with values between 0 and 2**niter-1)
     */
    std::vector<unsigned long long> coordinates_from_distance( unsigned long long h)
    {
        if (h > this->max_h)
            throw std::domain_error("h is greater than 2^(niter*dim)");
        auto x = this->_hilbert_integer_to_transpose(h);
        unsigned long long Z = 2ULL << (this->number_of_iterations -1);
        // Gray decode by H ^ (H/2)
        unsigned long long t = x[this->dimension-1] >> 1;
        for ( int i = this->dimension-1; i > 0; --i)
            x[i] ^= x[i-1];            
        x[0] ^= t;

        // Undo excess work
        unsigned long long Q = 2ULL;
        while (Q != Z)
        {
            unsigned long long P = Q - 1;
            for (int i = this->dimension-1; i >= 0; --i )
                if (x[i] & Q) // invert
                    x[0] ^= P;
                else
                { // Exchange
                    unsigned long long t = (x[0] ^ x[i])&P;
                    x[0] ^= t;
                    x[i] ^= t;
                }
            Q <<= 1;
        }
        // done
        return x;
    }
/*    def distance_from_coordinates(self, x_in):
        """Return the hilbert distance for a given set of coordinates.

        Args:
            x_in (list): transpose of h*/

    /**
     * @brief      Return the hilbert distance for a given set of coordinates.
     *
     * @param[in]  x_in  transpose of h
                         (n components with values between 0 and 2**p-1)
     *
     * @return     integer distance along hilbert curve
     */
    unsigned long long distance_from_coordinates( const std::vector<unsigned long long>& x_in)
    {
        if (x_in.size() != this->dimension)
            throw std::domain_error("x must have size equal to the dimension of space");
        for ( const auto& val : x_in )
            if (val > this->max_x)
                throw std::domain_error("A value in x_in has invalid coordinate");
        unsigned long long M = 1ULL << (this->number_of_iterations -1);

        std::vector<unsigned long long> x(x_in);
        // Inverse undo excess work
        unsigned long long Q = M;
        while (Q > 1ULL)
        {
            unsigned long long P = Q - 1;
            for ( int i = 0; i < this->dimension; ++i )
                if (x[i] & Q)
                    x[0] ^= P;
                else
                {
                    unsigned long long t = (x[0] ^ x[i]) & P;
                    x[0] ^= t;
                    x[i] ^= t;
                }
            Q >>= 1;
        }
        // Gray encode
        for ( int i = 1; i < this->dimension; ++i )
            x[i] ^= x[i-1];
        unsigned long long t = 0;
        Q = M;
        while (Q > 1)
        {
            if (x[this->dimension-1] & Q)
                t ^= Q - 1;
            Q >>= 1;
        }
        for (int i = 0; i < this->dimension; ++i )
            x[i] ^= t;

        unsigned long long h = this->_transpose_to_hilbert_integer(x);
        return h;
    }

private:
    /**
     * @brief      Convertit un entier de Hilbert en une coordonnée à n dimensions (valeurs entre 0 et 2^p-1)
     *
     * @param[in]  h     L'entier d'Hilbert
     *
     * @return     La coordonnée du point entier sur l'hypercube [0,2^p-1]^n
     */
    std::vector<unsigned long long> _hilbert_integer_to_transpose( unsigned long long h )
    {
        std::vector<unsigned long long> x;
        for ( unsigned idim = 0; idim < this->dimension; ++idim )
        {
            x[idim] = 0;
            for ( unsigned iter = 0; iter < this->number_of_iterations; ++iter )
            {
                if (h & (1ULL<<(idim+iter*this->dimension)))
                    x[idim] += (1ULL<<iter);
            }
        }
        //h_bit_str = _binary_repr(h, self.p*self.n)
        //x = [int(h_bit_str[i::self.n], 2) for i in range(self.n)]
        return x;
    }

    /**
     * @brief      Restaure un entier d'hilbert de sa transposée x
     *
     * @param[in]  x     Les coordonnées de point d'hilbert
     *
     * @return     L'entier d'hilbert (distance dans la courbe d'hilbert)
     */
    unsigned long long _transpose_to_hilbert_integer( const std::vector<unsigned long long>& x )
    {
        unsigned long long h = 0ULL;
        for ( unsigned idim = 0; idim < this->dimension; ++idim)
        {
            for ( unsigned iter = 0; iter < this->number_of_iterations ; ++iter )
            {
                if (x[idim] & (1ULL<<iter))
                {
                    h |= (1ULL<<((this->dimension-idim-1)+iter*this->dimension));
                    //std::cout << "1";
                }
                //else std::cout << "0";
            }
            //std::cout << " ";
        }
    /*def _transpose_to_hilbert_integer(self, x):
        x_bit_str = [_binary_repr(x[i], self.p) for i in range(self.n)]
        h = int(''.join([y[i] for i in range(self.p) for y in x_bit_str]), 2)
        return h
*/
        return h;
    }

    unsigned number_of_iterations;// Nbre iterations utilise dans la courbe d'Hilbert
    unsigned dimension;           // Dimension de l'espace
    unsigned long long max_h;     // Distance maximale possible pour cette courbe d'Hilbert
    unsigned long long max_x;     // Valeur max d'une coordonnee dans toute direction
};


#endif