/*=========================================================================
 
 Program:   itkSwathChainCodePathFilter.h
 Language:  C++
 Date:      $Date: 2010/04/13 15:16:04 $
 Version:   $Revision: 1.3 $
 
 Copyright (c) John Galeotti & George Stetten
 
 This software is distributed WITHOUT ANY WARRANTY; without even 
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 PURPOSE.
 
 =========================================================================*/

#ifndef __itkSwathChainCodePathFilter_h
#define __itkSwathChainCodePathFilter_h

#include "itkPathAndImageToPathFilter.h"
#include <list> // used by the private funciton FindPstepNeighborhoodOffsets
#include "itkConstantBoundaryCondition.h"

namespace itk
{
  
/** \class SwathChainCodePathFilter
 * \brief Filter that optimizes a path relative to an image in ND.
 *
 * SwathChainCodePathFilter produces a minimally-connected (vertex-connected)
 * chain code representation of a path that is optimal with respect to an image
 * and an original maximally-connected chain code path (sometimes referred to as
 * an "initial contour").  A hypercube "fovea" traces along the initial contour,
 * sweeping out a "swath" through the image. Dynamic programming is used to find
 * the "optimal" path through the part of the image swept out by the swath,
 * where "optimal" is defined as the path with the maximum sum of merits for its
 * individual steps.  A user-specified metrit function (not to be confused with
 * registration metrics) is responsible for evaluating the merit of an
   individual chaincode step located at a specific point in the input image:
 *
 * merit( direction of step, input-image index of end of step, input image)
 *
 * The optimal path is constrained to have the same number of steps as the
 * original path, but the "physical" length of the path may change since
 * diagonal non-maximally-connected steps have longer "physical" length than do
 * maximally-connected steps.
 *
 * \ingroup PathFilters
 */
template <class TChainCodePath, class TImage>
class ITK_EXPORT SwathChainCodePathFilter : public
PathAndImageToPathFilter< TChainCodePath, TImage, TChainCodePath >
{
public:
  /** Standard class typedefs. */
  typedef SwathChainCodePathFilter                            Self;
  typedef PathAndImageToPathFilter< TChainCodePath, TImage,
                                            TChainCodePath >  Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);
  
  /** Run-time type information (and related methods). */
  itkTypeMacro(SwathChainCodePathFilter, PathAndImageToPathFilter);

  /** ImageDimension constant */
  itkStaticConstMacro(ImageDimension, unsigned int, TImage::ImageDimension);
  
  /** Some convenient typedefs. */
  typedef          TChainCodePath           PathType;
  typedef typename PathType::Pointer        PathPointer;
  typedef typename PathType::ConstPointer   PathConstPointer;
  typedef typename PathType::InputType      PathInputType;
  
  typedef TImage                            ImageType;
  typedef typename ImageType::ConstPointer  ImageConstPointer;
  
  typedef typename PathType::IndexType      IndexType;
  typedef typename PathType::OffsetType     OffsetType;
  typedef typename ImageType::SizeType      SizeType;
  typedef typename ImageType::RegionType    RegionType;
  
  typedef typename std::list<OffsetType>    OffsetListType;
  typedef typename OffsetListType::iterator OffsetListIteratorType;
  
  
  
  /** Complex data structure typedefs. */
  
  // MeritImageImageType represents an image of pointers to images of doubles
  typedef Image< double,
    itkGetStaticConstMacro(ImageDimension) >      MeritImageType;
  typedef typename MeritImageType::Pointer        MeritImagePointer;
  typedef Image< MeritImagePointer,
    itkGetStaticConstMacro(ImageDimension) >      MeritImageImageType;
  typedef typename MeritImageImageType::Pointer   MeritImageImagePointer;
  
  // OffsetImageImageType represents an image of pointers to images of offsets
  typedef Image< OffsetType,
    itkGetStaticConstMacro(ImageDimension) >      OffsetImageType;
  typedef typename OffsetImageType::Pointer       OffsetImagePointer;
  typedef Image< OffsetImagePointer,
    itkGetStaticConstMacro(ImageDimension) >      OffsetImageImageType;
  typedef typename OffsetImageImageType::Pointer  OffsetImageImagePointer;
  
  typedef ConstantBoundaryCondition<MeritImageType> MeritBoundaryConditionType;
  
  
  /** Set the Fovea's radius.
   * Each side of the fovea hypercube will be of length 2*FoveaRadius+1 */
  itkSetMacro( FoveaRadius, unsigned int )
  
  /** Set whether or not to automatically trim the resultant output path.
   * By default, this is set to off (i.e. trim the path). */
  itkSetMacro( DoNotTrimPath, bool )
  itkBooleanMacro( DoNotTrimPath )
    
    
  /** Set the filter's merit-function pointer.  Invoke as
   * SetMeritFunctionPointer( MeritFunctionNameWithoutParentheses ) */
  virtual void SetMeritFunctionPointer( double (* MeritFunctionPointer)
    (const OffsetType step, const IndexType index, ImageConstPointer image) )
    {
    itkDebugMacro("setting MeritFunctionPointer to " << MeritFunctionPointer );
    if (this->m_MeritFunctionPointer != MeritFunctionPointer)
      {
      this->m_MeritFunctionPointer = MeritFunctionPointer;
      this->Modified();
      }
    }
  
protected:
  SwathChainCodePathFilter();
  virtual ~SwathChainCodePathFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData(void);

private:
  SwathChainCodePathFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  // each side of the fovea hypercube will be of length 2*m_FoveaRadius+1
  unsigned int m_FoveaRadius;
  
  bool m_DoNotTrimPath; // if set to true, then do not try to trim the path; false by default
  
  // pointer to the merit function
  // this should only be accessed through the MeritFunction() wrapper below
  double (* m_MeritFunctionPointer)(const OffsetType  step,
                                    const IndexType   index,
                                    ImageConstPointer image);
  
  // Wrapper function used to access m_MeritFunctionPointer
  virtual inline double MeritFunction(  const OffsetType  step,
    const IndexType   index, ImageConstPointer image)
    {
    // We really should do a lot of cacheing here...
    return (* m_MeritFunctionPointer)(step, index, image);
    }
  
  // Find the offsets to generate a shaped neighborhood coresponding to the
  // starting locations of all psteps coresponding to the current step of the
  // original path
  void FindPstepNeighborhoodOffsets( OffsetListType   &offsetList,
                                     const OffsetType  step );
  
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSwathChainCodePathFilter.txx"
#endif

#endif
