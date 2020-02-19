// A child wellBase class
// providing an implementation of pre-&-post correctors
class childWell : public base2IsoWell
{
public:
    TypeName("childWell")
    childWell(
            const word& name,
            const dictionary& wellDict,
            const fvMesh& mesh,
            const base2IsoWellModel& corrector
    ): base2IsoWell(name, wellDict, mesh, corrector){}
    virtual ~childWell(){}
    virtual bool writeData(Ostream&) const {}
    virtual void preCorrect(){ Info << "Executing preCorrect"  << endl; }
    virtual void postCorrect(){ Info << "Executing postCorrect"  << endl; }
};
defineTypeNameAndDebug(childWell, 0);
addToRunTimeSelectionTable(base2IsoWell, childWell, dictionary);
