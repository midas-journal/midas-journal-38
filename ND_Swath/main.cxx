#include <iostream.h>
#include <exception>
#include "itkExceptionObject.h" // only needed for extra ITK output

extern int itkSwathChainCodePathFilterTest(int, char**);


int main(int argc, char *argv[])
{
  int status=0;
  
try
  {
  status = status || itkSwathChainCodePathFilterTest(argc, argv);
  cout << "\n----------\n" << endl;
  }
catch ( itk::ExceptionObject& e )
  {
  std::cerr << e;
  std::cerr << "Done printing itk::ExceptionObjects" << std::endl;
  status = 1;
  }
catch ( std::exception& e )
  {
  std::cerr << "Exception: " << e.what();
  std::cerr << "Done printing std::Exceptions" << std::endl;
  status = 1;
  }
  cout << "Status = " << status << endl;
  return status;
}
