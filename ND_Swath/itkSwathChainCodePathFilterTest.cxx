/*=========================================================================

  Program:   itkSwathChainCodePathFilterTest.cxx
  Language:  C++
  Date:      $Date: 2010/04/13 15:16:04 $
  Version:   $Revision: 1.4 $

  Copyright (c) John Galeotti & George Stetten

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkPolyLineParametricPath.h"
#include "itkChainCodePath.h"
#include "itkFourierSeriesPath.h"
#include "itkPathToChainCodePathFilter.h"
#include "itkChainCodeToFourierSeriesPathFilter.h"

#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkDiscreteGaussianImageFilter.h"

#include "itkSwathChainCodePathFilter.h"

#include "itkPathToImageFilter.h"
#include "itkImageFileWriter.h"


double meritFunction(const itk::Offset<2> step, const itk::Index<2> index,
                     itk::Image<double,2>::ConstPointer image);


int itkSwathChainCodePathFilterTest(int, char*[])
{
  typedef itk::Image<unsigned char, 2>                ImageType;         
  typedef itk::Image<double, 2>                       FloatImageType;
  typedef itk::PolyLineParametricPath<2>              InPathType;        
  typedef itk::ChainCodePath<2>                       ChainPathType;     
  typedef itk::FourierSeriesPath<2>                   FSPathType;        

  typedef InPathType::VertexType                      VertexType;        
  typedef InPathType::OffsetType                      OffsetType;        
  typedef InPathType::InputType                       InPathInputType;   

  typedef ImageType::IndexType                        IndexType;         
  typedef ImageType::SizeType                         SizeType;          

  typedef ChainPathType                               OutputPathType;    

  // pre-process the image
  typedef itk::RescaleIntensityImageFilter<ImageType,FloatImageType>      CastFilterType;
  typedef itk::DiscreteGaussianImageFilter<FloatImageType,FloatImageType> SmoothFilterType;
  
  // pre-process the path
  typedef itk::PathToChainCodePathFilter<InPathType,ChainPathType>          PathFilter1Type;
  typedef itk::ChainCodeToFourierSeriesPathFilter<ChainPathType,FSPathType> PathFilter2Type;
  
  // test the filter
  typedef itk::SwathChainCodePathFilter <OutputPathType,FloatImageType>     TestFilterType;
  
  // save the results
  typedef itk::PathToImageFilter< ChainPathType, ImageType >          Output1FilterType;
  typedef itk::PathToImageFilter< OutputPathType, ImageType >         Output2FilterType;
  typedef itk::RescaleIntensityImageFilter<FloatImageType,ImageType>  Output3FilterType;
  
  
  std::cerr<<"SwathChainCodePathFilter"<<std::endl;
  
  
  
  // Setup the path
  std::cout << "Making a square Path with v0 at (17,17) -> (13,51) -> (50,50) -> (49,13)" << std::endl;
  InPathType::Pointer inPath;
  VertexType        v;
  inPath = InPathType::New();
  v.Fill(17);
  inPath->AddVertex(v);
  v[0]=13;
  v[1]=51;
  inPath->AddVertex(v);
  v.Fill(50);
  inPath->AddVertex(v);
  v[0]=49;
  v[1]=13;
  inPath->AddVertex(v);
  v.Fill(17);
  inPath->AddVertex(v);
  
  // Setup the first path filter
  PathFilter1Type::Pointer pathFilter1 = PathFilter1Type::New();
  pathFilter1->SetInput(inPath);
  
  // Setup the second path filter
  PathFilter2Type::Pointer pathFilter2 = PathFilter2Type::New();
  pathFilter2->SetInput(pathFilter1->GetOutput());
  pathFilter2->SetNumberOfHarmonics(7); // make a nice, round, path for the swath
  
  
  
  // Setup the image
  std::cout << "Making a 32x32 black square centered in a 64x64 white image"<<std::endl;
  // a trace-around path should go from 15,15 -> 48,15 -> 48,48 -> 15,48 -> 15,15
  ImageType::Pointer  inImage = ImageType::New();
  IndexType start;
  start[0]=0;
  start[1]=0;
  ImageType::SizeType size;
  size[0]=64;
  size[1]=64;
  ImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);
  inImage->SetRegions(region);
  double spacing[ ImageType::ImageDimension ];
  spacing[0]=1.0;
  spacing[1]=1.0;
  inImage->SetSpacing(spacing);
  inImage->Allocate();
  typedef itk::ImageRegionIterator<ImageType> ImageItType;
  ImageItType it( inImage, inImage->GetRequestedRegion() );
  it.GoToBegin();
  IndexType pixelIndex;
  while( !it.IsAtEnd() )
    {
    pixelIndex = it.GetIndex();
    if( pixelIndex[0] >= int(size[0]/4) && pixelIndex[0] < int(size[0]*3/4) &&
        pixelIndex[1] >= int(size[1]/4) && pixelIndex[1] < int(size[1]*3/4) )
      it.Set(0);
    else
      it.Set(255);
    ++it;
    }
  
  // Cast the input image into a double image
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput( inImage );
  castFilter->SetOutputMinimum(0);
  castFilter->SetOutputMaximum(1.0);

  // Smooth the (double pixel type) input image
  SmoothFilterType::Pointer smoothFilter = SmoothFilterType::New();
  smoothFilter->SetInput( castFilter->GetOutput() );
  double gaussianVariance = 1.0;
  // We want a fast 3x3 kernel.  Gausian Operator will not truncate its kernel
  // width to any less than a 5x5 kernel (kernel width of 3 for 1 center pixel +
  // 2 edge pixels).  However, GausianOperator always uses at least a 3x3
  // kernel, and so setting the maximum error to 1.0 (no limit) will make it
  // stop growing the kernel at the desired 3x3 size.
  double maxError = 0.9;  
  smoothFilter->SetUseImageSpacingOff();
  smoothFilter->SetVariance( gaussianVariance );
  smoothFilter->SetMaximumError( maxError );
  
  
  
  // Setup the test filter
  std::cerr << "Creating the test filter" << std::endl;
  TestFilterType::Pointer testFilter = TestFilterType::New();
  std::cerr << "Setting up the test filter" << std::endl;
  testFilter->SetPathInput(pathFilter1->GetOutput());
  //testFilter->SetImageInput(castFilter->GetOutput());
  testFilter->SetImageInput(smoothFilter->GetOutput());
  OutputPathType::Pointer outPath = testFilter->GetOutput();
  testFilter->SetFoveaRadius(4);  // a 7x7 fovea
  testFilter->SetMeritFunctionPointer( meritFunction );
  
  // Setup the input & output path images
  std::cerr << "Setting up the output path images" << std::endl;
  Output1FilterType::Pointer output1Filter = Output1FilterType::New();
  Output2FilterType::Pointer output2Filter = Output2FilterType::New();
  size[0]=64;
  size[1]=64;
  output1Filter->SetSize(size);  // same size as the input image
  output2Filter->SetSize(size);  // same size as the input image
  output1Filter->SetPathValue(255);
  output2Filter->SetPathValue(255);
  output1Filter->SetInput( testFilter->GetPathInput() );
  output2Filter->SetInput( testFilter->GetOutput() );
  ImageType::Pointer inPathImage = output1Filter->GetOutput();
  ImageType::Pointer outImage =    output2Filter->GetOutput();
  
  // Setup the test-filter-input-image for writting
  Output3FilterType::Pointer output3Filter = Output3FilterType::New();
  //output3Filter->SetInput( smoothFilter->GetOutput() );
  output3Filter->SetInput( testFilter->GetImageInput() );
  output3Filter->SetOutputMinimum(0);
  output3Filter->SetOutputMaximum(255);
  ImageType::Pointer swathInputImage = output3Filter->GetOutput();
  
  
  // Testing PrintSelf
  std::cerr << "testFilter is:\n" << testFilter << std::endl;
  
  // Update the pipeline
  std::cerr << "Running the Pipeline: ";
  //testFilter->DebugOn();
  outImage->Update();
  //testFilter->DebugOff();
  std::cerr << "[DONE]" << std::endl;
  
  // Save the output images
  std::cerr << "Saving images"<<std::endl;
  itk::ImageFileWriter<ImageType>::Pointer writer
    = itk::ImageFileWriter<ImageType>::New();
  writer->SetFileName( "SwathFilterPathOut.png" );
  output2Filter->Update();
  writer->SetInput( output2Filter->GetOutput() );
  writer->Write();
  writer->SetFileName( "SwathFilterPathIn.png" );
  output1Filter->Update();
  writer->SetInput( output1Filter->GetOutput() );
  writer->Write();
  writer->SetFileName( "SwathFilterInputImage.png" );
  output3Filter->Update();
  writer->SetInput( output3Filter->GetOutput() );
  writer->Write();
  
  for( unsigned int i=0; i<outPath->NumberOfSteps(); i++ )
  std::cout<<"OptimalStep["<<i<<"] = "<<outPath->Evaluate(i)<<std::endl;
  
  
  // Now find the output image w/o path trimming
  std::cerr << "\n----------\n";
  std::cerr << "Turning off path trimming." << std::endl;
  testFilter->DoNotTrimPathOn();
  // Update the pipeline
  std::cerr << "Running the Pipeline: ";
  //testFilter->DebugOn();
  outImage->Update();
  //testFilter->DebugOff();
  std::cerr << "[DONE]" << std::endl;
  // Save the output image
  std::cerr << "Saving image"<<std::endl;
  writer->SetFileName( "SwathFilterRawPathOut.png" );
  output2Filter->Update();
  writer->SetInput( output2Filter->GetOutput() );
  writer->Write();

  
  return EXIT_SUCCESS;
}



// Setup a merit function to trace the exterior bounding pixels of a "blob"
#include "vnl/vnl_math.h"
double meritFunction(const itk::Offset<2> step, const itk::Index<2> index,
                     itk::Image<double,2>::ConstPointer image)
{
  const double interiorValue=0.0;  // desired pixel value for interior pixels
  const double exteriorValue=1.0;  // desired pixel value for exterior pixels
  const bool clockwise = true;  // does the path trace CW or CCW around the blob?
  
  // static lookup tables
  static int freemanCode[3][3]; // freemanCode[ offset[0] + 1 ][ offset[1] + 1 ]
  static itk::Offset<2> *reverseFreemanCode = NULL;  // stored in clockwise order
  
  
  // Freeman code representation of step
  int freemanStep;
  
  itk::Offset<2> boundaryOffset; // offset from the inside index to the outside index
  itk::Index<2>  insideIndex;    // index of the inside pixel
  double    insideValue;   // measured value of pixel that should be inside the blob
  double    outsideValue1; // measured value of pixel that should be outside the blob
  double    outsideValue2; // like outsideValue1, but for a second "outside" pixel
  
  
  
  // Initialize static lookup tables the first time this function is called
  if( !reverseFreemanCode )
    {
    reverseFreemanCode = new itk::Offset<2>[8];
    reverseFreemanCode[0][0] =  0;    reverseFreemanCode[0][1]=  1; // up
    reverseFreemanCode[1][0] =  1;    reverseFreemanCode[1][1]=  1;
    reverseFreemanCode[2][0] =  1;    reverseFreemanCode[2][1]=  0; // right
    reverseFreemanCode[3][0] =  1;    reverseFreemanCode[3][1]= -1;
    reverseFreemanCode[4][0] =  0;    reverseFreemanCode[4][1]= -1; // down
    reverseFreemanCode[5][0] = -1;    reverseFreemanCode[5][1]= -1;
    reverseFreemanCode[6][0] = -1;    reverseFreemanCode[6][1]=  0; // left
    reverseFreemanCode[7][0] = -1;    reverseFreemanCode[7][1]=  1;
    
    freemanCode[ reverseFreemanCode[0][0] + 1 ][ reverseFreemanCode[0][1] + 1 ] = 0;
    freemanCode[ reverseFreemanCode[1][0] + 1 ][ reverseFreemanCode[1][1] + 1 ] = 1;
    freemanCode[ reverseFreemanCode[2][0] + 1 ][ reverseFreemanCode[2][1] + 1 ] = 2;
    freemanCode[ reverseFreemanCode[3][0] + 1 ][ reverseFreemanCode[3][1] + 1 ] = 3;
    freemanCode[ reverseFreemanCode[4][0] + 1 ][ reverseFreemanCode[4][1] + 1 ] = 4;
    freemanCode[ reverseFreemanCode[5][0] + 1 ][ reverseFreemanCode[5][1] + 1 ] = 5;
    freemanCode[ reverseFreemanCode[6][0] + 1 ][ reverseFreemanCode[6][1] + 1 ] = 6;
    freemanCode[ reverseFreemanCode[7][0] + 1 ][ reverseFreemanCode[7][1] + 1 ] = 7;
    }
  
  // get the Freeman code representation of step
  freemanStep = freemanCode[ step[0]+1 ][ step[1]+1 ];

  // Are we dealing with a straight or diagonal step?
  if( 0==freemanStep%2 )
    {
    
    // step is straight:  check just one outside pixel value
    //------------------
    
    // get step rotated 90 deg. opposite the direction of the path around the blob
    boundaryOffset = ( clockwise ? reverseFreemanCode[(freemanStep+6)%8]
                                 : reverseFreemanCode[(freemanStep+2)%8] );
    insideIndex   = index - boundaryOffset;
    
    insideValue   = image->GetPixel( insideIndex );
    outsideValue1 = image->GetPixel( index );
    
    // penalize for differences in both inside and outside values
    return -( vnl_math_abs( insideValue   - interiorValue ) +
              vnl_math_abs( outsideValue1 - exteriorValue ) );
    }
  else
    {
    
    // step is diagonal:  check two outside pixel values
    //------------------
    
    // get step rotated 45 deg. opposite the direction of the path around the blob
    boundaryOffset = ( clockwise ? reverseFreemanCode[(freemanStep+7)%8]
                                 : reverseFreemanCode[(freemanStep+1)%8] );
    insideIndex   = index - boundaryOffset;
    // get step rotated 90 deg. opposite the direction of the path around the blob
    boundaryOffset = ( clockwise ? reverseFreemanCode[(freemanStep+6)%8]
                                 : reverseFreemanCode[(freemanStep+2)%8] );
    
    insideValue   = image->GetPixel( insideIndex );
    outsideValue1 = image->GetPixel( index );
    outsideValue2 = image->GetPixel( insideIndex + boundaryOffset );
    // penalize for differences in both inside and outside values
    // *** Is 0.5 the proper multiplier, or is such scaling appropriate? ***
    return -(     vnl_math_abs( insideValue   - interiorValue ) +
              0.5*vnl_math_abs( outsideValue1 - exteriorValue ) +
              0.5*vnl_math_abs( outsideValue2 - exteriorValue ) );
    }
  
}
