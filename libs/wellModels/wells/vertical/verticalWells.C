
#include "vertical.C"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace wells
{
    // Isotropic and anisotropic spesializations 
    typedef vertical<volScalarField> isoVertical;
    typedef vertical<volTensorField> anisoVertical;

    // Member functions specialization for isotropic media

        //- Estimate well equivalent radius
        template<>
        void isoVertical::estimateEquivRadius(const Foam::volScalarField& k)
        {
            // Isotropic media ---> Re depends only on geometry
            BHPData_.re.setSize(wellSet_.toc().size());
            forAll(wellSet_, celli)
            {
                BHPData_.re[celli] = 
                    (0.14/0.5)
                    *sqrt(pow(BHPData_.h[celli][0],2)+pow(BHPData_.h[celli][1],2))
                    /2.0;
            }
        }

        //- Return sqrt(K11*K22)
        template<>
        scalar isoVertical::sqrtKK(const Foam::volScalarField& k, const label& cellID)
        {
            return k[cellID];
        }

        //- Return coefficients for injection rates of multi-cell wells
        template<>
        scalarList isoVertical::correctInjectionIntoCells(const Foam::volScalarField& k)
        {
            if (not(weightCellRates_)) return List<scalar>();
            scalarList theCoeffs(wellSet_.toc().size(),1.0/wellSet_.toc().size());
            // Read from dictionary only if specified
            // Otherwise return theCoeffs
            return wellDict_.lookupOrDefault<List<scalar> >("cellCoeffs",theCoeffs);
        }

    // Member functions specialization for anisotropic media

        //- Estimate well equivalent radius
        template<>
        void anisoVertical::estimateEquivRadius(const Foam::volTensorField& k)
        {
            // Anisotropic media 
            BHPData_.re.setSize(wellSet_.toc().size());
            forAll(wellSet_, celli)
            {
                label cellID = wellSet_.toc()[celli];
                scalar k11 = k[cellID].xx();
                scalar k22 = k[cellID].yy();
                BHPData_.re[celli] = 
                    (0.14/0.5)
                    *sqrt(sqrt(k22/k11)*pow(BHPData_.h[celli][0],2)+sqrt(k11/k22)*pow(BHPData_.h[celli][1],2))
                    /(pow(k11/k22,0.25)+pow(k22/k11,0.25));
            }
        }

        //- Return sqrt(K11*K22)
        template<>
        scalar anisoVertical::sqrtKK(const Foam::volTensorField& k, const label& cellID)
        {
            return sqrt(k[cellID].xx()*k[cellID].yy());
        }

        //- Return coefficients for injection rates of multi-cell wells
        template<>
        scalarList anisoVertical::correctInjectionIntoCells(const Foam::volTensorField& k)
        {
            // Weight Injection volumes per cell with (K.V)/tV
            // This should be efficient because the number of cells is small
            scalarList theCoeffs(wellSet_.toc().size(),0.0);
            forAll(theCoeffs,celli)
            {
                scalar kc = k[celli][0]*k[celli][4]; // Initialize to Kxx.Kyy
                if (k[celli][0]< VSMALL) kc = k[celli][4]*k[celli][4]; 
                if (k[celli][4]< VSMALL) kc = k[celli][0]*k[celli][0]; 
                theCoeffs[celli] = sqrt(kc)*mesh_.V()[celli]/tV_;
            }
            scalar cCoeffs = sum(theCoeffs);
            forAll(theCoeffs,celli)
            {
                theCoeffs[celli] /= (cCoeffs+VSMALL);
            }

            return theCoeffs;
        }

    // Add both versions to RTS table
    defineTemplateTypeNameAndDebug(isoVertical, 0);
    defineTemplateTypeNameAndDebug(anisoVertical, 0);
    
    addToRunTimeSelectionTable
    (
        well,
        isoVertical,
        dictionary
    );
    addToRunTimeSelectionTable
    (
        well,
        anisoVertical,
        dictionary
    );
}
}
