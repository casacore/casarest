#ifdef HAVE_ADIOS2
#include <synthesis/Images/ADIOSImage.h>
#endif


int main(int argc, char** argv)
{
#ifdef HAVE_ADIOS2
    casacore::ADIOSImage<float> adios("test");
#endif
    return 0;
}
