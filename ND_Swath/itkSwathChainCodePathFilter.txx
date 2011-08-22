/*=========================================================================

  Program:   itkSwathChainCodePathFilter.txx
  Language:  C++
  Date:      $Date: 2010/04/13 15:16:04 $
  Version:   $Revision: 1.4 $

  Copyright (c) John Galeotti & George Stetten

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.

=========================================================================*/

#ifndef _itkSwathChainCodePathFilter_txx
#define _itkSwathChainCodePathFilter_txx

#include "itkImageRegionIterator.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkConstNeighborhoodIterator.h" // follows the path
#include "itkConstShapedNeighborhoodIterator.h"

#include "itkNumericTraits.h"

#include "itkSwathChainCodePathFilter.h"


namespace itk
{

/**
 * Constructor
 */
template <class TChainCodePath, class TImage>
SwathChainCodePathFilter<TChainCodePath, TImage>
::SwathChainCodePathFilter()
{
  m_FoveaRadius = 2;
  m_DoNotTrimPath = false;
  m_MeritFunctionPointer = NULL;
}


/**
 * PrintSelf
 */
template <class TChainCodePath, class TImage>
void
SwathChainCodePathFilter<TChainCodePath, TImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "FoveaRadius:  " << m_FoveaRadius << std::endl;
  os << indent << "DoNotTrimPath:  " << m_DoNotTrimPath << std::endl;
  os << indent << "MeritFunctionPointer:  " << m_MeritFunctionPointer << std::endl;
}



/**
 * GenerateData Performs the reflection
 */
template <class TChainCodePath, class TImage>
void
SwathChainCodePathFilter<TChainCodePath, TImage>
::GenerateData( void )
{
  // ----------------------------------------------------------------------- //
  //   Setup basic data types
  // ----------------------------------------------------------------------- //
  itkDebugMacro( << "GenerateData():  Initializing simple data structures." );
  
  ImageConstPointer inputImage  = this->GetImageInput();
  PathConstPointer  inputPath   = this->GetPathInput();
  PathPointer       outputPath  = this->GetOutput();
  unsigned int      pathLength  = inputPath->NumberOfSteps();
  SizeType          foveaRadius; foveaRadius.Fill(m_FoveaRadius);
  SizeType          foveaSize;   foveaSize.Fill(2*m_FoveaRadius+1);
  RegionType        foveaRegion;
    {
    IndexType start;
    start.Fill(0);
    foveaRegion.SetIndex(start);
    foveaRegion.SetSize(foveaSize);
    }
  OffsetListType    dirtyPath;  // dirtyPath is cleaned up then copied to outputPath.
  bool              openPath;


  
  
  // ----------------------------------------------------------------------- //
  //   Create and initialize the very complex data structures.
  // ----------------------------------------------------------------------- //
  
  // Some Notes on Notation:
  //   
  //   A path is an entire, complete chain code
  //   
  //   A step is part of a path and is stored as an offset.
  //   
  //   A trail is an initial segment of a path (which may reach all the way to
  //     the end of the path).  Trails are more of a concept and are not stored
  //     as chain codes, but a trail can be converted to a chain code by
  //     iterating through some data structures. There exists one optimal trail
  //     from each possible starting point of a path to each possible ending
  //     point of a path.  By growing a path from a single step to its full
  //     length, the coresponding set of optimal trails can also be grown until
  //     every optimal possible path from every possible starting point to every
  //     possible ending point is represented by a trail.  The best trail from
  //     this set of optimal trails can then be converted to a path from the
  //     optimal starting location to the optimal ending location following the
  //     optimal course in between.
  //   
  //   A tstep (trail step) is a part of a trail and is stored as an offset.  A
  //     trail is converted to a path by concatenating its tsteps (which are not
  //     generally stored in order due to the nature of dynamic programming).
  //   
  //   A pstep (possible step) has the potential to be a tstep (and hence the
  //     possiblility of eventually becomming a step of an actual path).  For a
  //     given trail starting location and a given ending location, each valid
  //     pstep "pointing" to that ending location can be evaluated for both its
  //     own merit as well as for the merit of the trail beginning at the trail
  //     starting location and ending where the pstep begins (the trail that the
  //     pstep may be concatenated to if it were chosen to become a tstep). 
  //     Psteps are also stored as offsets. 
  // 
  
  // MeritImage and tstepImage each have one index for every trail from a
  // single starting location to every possible ending location.  We must keep
  // track of every starting location, however, and so we must maintain a
  // meritImage and a tstepImage for every trail starting location.  The number
  // of trails from a single starting location = the number of starting
  // locations = the number of indices in the fovea.  Therefore, the number of
  // trails we must keep track of is equal to the square of the number of
  // indices in the fovea and so meritImage, meritImageSet, tstepImage, and
  // tstepImageImage are each of the same size as the fovea.
  // 
  // As the trails grow to the length of the input path, we must keep track of
  // the continually added merits and offsets.  Because we use dynamic
  // programming, we only need to keep track of the merits for the current
  // length trails and for the previous (1 less) length trails.  However, to be
  // able to reconstruct the optimal path from the optimal trail, we must keep
  // track of all the tsteps added to each trail for the entire length of every
  // trail (which will eventually grow to the length of the input path).
  // 
  // **************************************************************************
  // Therefore, the 3 top level data structures are:
  // 
  //   meritImageSet & previousMeritImageSet:
  //     images of smartpointers to images
  //     store the merits of all trails of current length and preceeding length
  //   
  //   tstepImageSet[]:
  //     array of images of smartpointers to images
  //     store the final tsteps of all optimal trails for every length of trail
  // **************************************************************************
  //
  MeritImagePointer       meritImage;       // image of trail merits (doubles)
  MeritImageImagePointer  meritImageSet;    // image of pointers to meritImages
  MeritImagePointer       previousMeritImage;     // just like meritImage
  MeritImageImagePointer  previousMeritImageSet;  // just like meritImageSet
  OffsetImagePointer      tstepImage;       // image of tsteps (offsets)
  OffsetImageImagePointer tstepImageImage;  // image of pointers to tstepImages
  OffsetImageImagePointer *tstepImageSet;   // array of tstepImageImage; set of
                                            // all tsteps for each position
                                            // along each trail
  
  // Initialize the image-based data structures
    {
    itkDebugMacro( << "GenerateData():  Initializing complex data structures." );
    
    // Create the data structures of merits
      {
      meritImageSet = MeritImageImageType::New();
      meritImageSet->SetRegions(foveaRegion);
      meritImageSet->Allocate();
      previousMeritImageSet = MeritImageImageType::New();
      previousMeritImageSet->SetRegions(foveaRegion);
      previousMeritImageSet->Allocate();

      typedef ImageRegionIteratorWithIndex<MeritImageImageType> iteratorType;
      // fill meritImageSet
        {
        iteratorType it( meritImageSet, foveaRegion );
        for( it.GoToBegin(); !it.IsAtEnd(); ++it )
          {
          meritImage = MeritImageType::New();
          meritImage->SetRegions(foveaRegion);
          meritImage->Allocate();
          it.Set(meritImage);
          }
        }
      // fill previousMeritImageSet with INITIALIZED meritImages
        {
        IndexType startingIndex;
        // iterate over all starting indices for trails
        iteratorType it( previousMeritImageSet, foveaRegion );
        for( it.GoToBegin(); !it.IsAtEnd(); ++it )
          {
          meritImage = MeritImageType::New();
          meritImage->SetRegions(foveaRegion);
          meritImage->Allocate();
          
          // Initialize all merit values so that each "zero-length previous"
          // trail has the same starting and ending location.
          startingIndex = it.GetIndex();
          // iterate over all ending indices for 0-length trails
          ImageRegionIteratorWithIndex<MeritImageType> it2(meritImage,foveaRegion);
          for( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 )
            {
            if( it2.GetIndex() == startingIndex )
              {
              // this represents a valid trail configuration
              it2.Set(NumericTraits<double>::Zero);
              }
            else
              {
              // this invalid 0-length trail must never be used:
              it2.Set(-NumericTraits<double>::infinity());
              }
            }
          
          it.Set(meritImage);
          }
        }
      }

    // Create the data structure of offsets
      {
      tstepImageSet = new OffsetImageImagePointer[ pathLength ];
      
      for( unsigned int distanceDownTrails=0; distanceDownTrails < pathLength;
           distanceDownTrails++ )
        {
        tstepImageImage = OffsetImageImageType::New();
        tstepImageImage->SetRegions(foveaRegion);
        tstepImageImage->Allocate();

        ImageRegionIterator<OffsetImageImageType> it( tstepImageImage, foveaRegion );
        for( it.GoToBegin(); !it.IsAtEnd(); ++it )
          {
          tstepImage = OffsetImageType::New();
          tstepImage->SetRegions(foveaRegion);
          tstepImage->Allocate();
          it.Set(tstepImage);
          }
        
        tstepImageSet[distanceDownTrails] = tstepImageImage;
        }
      }
    
    }
  
  
  
  // ----------------------------------------------------------------------- //
  //   Grow the optimal trails for all combinations of start and end indices
  // ----------------------------------------------------------------------- //
    {
    itkDebugMacro( << "GenerateData():  Growing trails." );
    
    OffsetType      step;
    OffsetListType  pstepNeighborhoodOffsetList;
    
    // create a neighborhood of trail-ending points called the fovea and center
    // it at the beginning of inputPath (at the beginning of the first step)
    ConstNeighborhoodIterator< ImageType > pathIt( foveaRadius, inputImage,
                                           inputImage->GetRequestedRegion() );
    pathIt.SetLocation( inputPath->GetStart() );  // fovea @ start of first step
    
    // iterate the fovea along the original path to grow the trails
    for( PathInputType pathInput = inputPath->StartOfInput();
         pathInput < inputPath->NumberOfSteps();
         ++pathInput )
      {
      itkDebugMacro(<<"GenerateData():  Growing length "<<pathInput+1<<" trails");
      
      step = inputPath->Evaluate( pathInput );    // get the current step
      pathIt += step;                             // fovea @ end of current step
      
      // What psteps could replace the current step?  Find the neighborhood
      // offset list for repeated use later.
      FindPstepNeighborhoodOffsets( pstepNeighborhoodOffsetList, step );
      
      // grow trails from every possible starting location
      ImageRegionIterator<MeritImageImageType>
        meritImageSetIterator( meritImageSet, foveaRegion );
      
      ImageRegionIterator<MeritImageImageType>
        previousMeritImageSetIterator( previousMeritImageSet, foveaRegion );
      
      ImageRegionIterator<OffsetImageImageType>
        tstepImageSetIterator( tstepImageSet[pathInput], foveaRegion );
      
      meritImageSetIterator.GoToBegin();
      previousMeritImageSetIterator.GoToBegin();
      tstepImageSetIterator.GoToBegin();
      
      while( !meritImageSetIterator.IsAtEnd() &&
             !previousMeritImageSetIterator.IsAtEnd() &&
             !tstepImageSetIterator.IsAtEnd() )
        {
        meritImage          = meritImageSetIterator.Get();
        previousMeritImage  = previousMeritImageSetIterator.Get();
        tstepImage          = tstepImageSetIterator.Get();
        
        //itkDebugMacro( << "GenerateData():  Growing length " << pathInput <<
        //               " trails from " << meritImageSetIterator.GetIndex() );

        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // Ok, at this point meritImage, previousMeritImage, & tstepImage are
        // all syncronized to a particular possible-trail-starting-location for
        // a particular step along the original path.
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        //
        
        
        
        // Prepare to iterate through every possible trail-ending-location.
        // Create tstepIt, meritIt, and previousMeritIt; pathIt already exists.
        
        ImageRegionIterator<OffsetImageType> tstepIt( tstepImage, foveaRegion );
        ImageRegionIterator<MeritImageType>  meritIt( meritImage, foveaRegion );
        
        // We must consider multiple psteps for each trail-ending-location; the
        // psteps are accessed by a neighborhood of starting-indices for psteps
        // all pointing to the current trail-ending-location.
        typename ConstShapedNeighborhoodIterator<MeritImageType>::RadiusType
          radiusOfOne; radiusOfOne.Fill(1);

        ConstShapedNeighborhoodIterator<MeritImageType>
          previousMeritIt( radiusOfOne, previousMeritImage, foveaRegion );
        
        // Make a neighborhood centered at the end of the previous tstep with a
        // pixel for each beginning index of each valid pstep.  The offsets to
        // generate the neighborhood were all precomputed for this step of the
        // original input path.  -Infinity must be returned if outside the
        // boundaries of previousMeritImage.
        previousMeritIt.ClearActiveList();
        for( OffsetListIteratorType it = pstepNeighborhoodOffsetList.begin();
          it != pstepNeighborhoodOffsetList.end();   it++ )
          {
          previousMeritIt.ActivateOffset( *it );
          }
        MeritBoundaryConditionType meritBoundaryCondition;
        meritBoundaryCondition.SetConstant( -NumericTraits<double>::infinity());
        previousMeritIt.OverrideBoundaryCondition( &meritBoundaryCondition );
        
        
        
        // Iterate through every possible trail-ending-location
        
        tstepIt.GoToBegin();
        meritIt.GoToBegin();
        previousMeritIt.GoToBegin();
        unsigned int pathItNeighborhoodIndex = 0;
        
        // all these ending conditions should be met at the same time
        //while( !tstepIt.IsAtEnd() &&
        //       !meritIt.IsAtEnd() &&
        //       !previousMeritIt.IsAtEnd() &&
        //       pathItNeighborhoodIndex < pathIt.Size() )
        while( pathItNeighborhoodIndex < pathIt.Size() )
          {
          //itkDebugMacro( << "GenerateData():  Growing length " << pathInput
          //               << " trail from " << meritImageSetIterator.GetIndex()
          //               << " to " << meritIt.GetIndex() );
          
          OffsetType  bestPstep = step;
          double      bestMerit = -NumericTraits<double>::infinity();
          double      oldMerit;
          double      newMerit;
          OffsetType  newPstep;
          
          // Find the best pstep to complete the current trail.
          // The best pstep will have the highest combination of its own merit
          // and the merit of the length n-1 trail that feeds it.
          // pstepIt iterates through the neighborhood of pstep tail positions.
          typename ConstShapedNeighborhoodIterator<MeritImageType>::
            ConstIterator pstepIt;
          for( pstepIt = previousMeritIt.Begin(); !pstepIt.IsAtEnd(); pstepIt++)
            {
            // if the feeding trail has -Infinity merit, don't bother
            oldMerit=pstepIt.Get();
            if( oldMerit != -NumericTraits<double>::infinity() )
              {
              // previousMeritIt iterates over the STARTING locations of the
              // psteps, so to find the actual offsets of the psteps represented
              // by the offsets of previousMeritIt, we need to subtract the
              // offsets of previousMeritIt from the offset of the original
              // step.
              newPstep = step - pstepIt.GetNeighborhoodOffset();
              newMerit = oldMerit + MeritFunction( newPstep,
                         pathIt.GetIndex(pathItNeighborhoodIndex), inputImage );
              
              // If this is the best merit we have seen so far, remember it;
              if( newMerit > bestMerit )
                {
                bestMerit = newMerit;
                bestPstep = newPstep;
                }
              }
            }
          
          tstepIt.Set( bestPstep );
          meritIt.Set( bestMerit );
        
          ++tstepIt;
          ++meritIt;
          ++previousMeritIt;
          ++pathItNeighborhoodIndex;
          }
        
        
        
        //
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        // done with this possible-trail-starting-location for this path step
        // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        
        ++meritImageSetIterator;
        ++previousMeritImageSetIterator;
        ++tstepImageSetIterator;
        }
      
      // We are done processing this step of the original path, so we need to
      // replace previousMeritImageSet with meritImageSet.  Because all the
      // values in meritImageSet will be overwritten on the next pass through
      // this loop, we can save the work of re-creating meritImageSet by just
      // swapping meritImageSet with previousMeritImageSet.
      MeritImageImagePointer tempMeritImageSet = previousMeritImageSet;
      previousMeritImageSet = meritImageSet;
      meritImageSet = tempMeritImageSet;
      
      } // done iterating along the path
      
      // Undo the merit image swapping done at the end of the last path step
      meritImageSet = previousMeritImageSet;
    }
  
  
  
  // ----------------------------------------------------------------------- //
  //   Find the optimal trail
  // ----------------------------------------------------------------------- //
    {
    itkDebugMacro( << "GenerateData():  Finding optimal trail." );
    
    IndexType bestTrailCurrentFoveaIndex;
    IndexType bestTrailStartFoveaIndex;
    IndexType bestTrailEndFoveaIndex;
    double trailMerit;
    double bestTrailMerit = -NumericTraits<double>::infinity();
    
    
    // Find the optimal trail
    
    // If the input path was closed, the output path must be closed too.
    openPath  = ( inputPath->EvaluateToIndex(inputPath->EndOfInput())  !=
                  inputPath->EvaluateToIndex(inputPath->StartOfInput())   );
    itkDebugMacro( << "GenerateData():  Finding optimal trail: Is path open?: "<<openPath );
    
    // Iterate over all trails (from all starting to all ending locations)
    // looking for the best trail
    ImageRegionIteratorWithIndex<MeritImageImageType>
      meritImageSetIterator( meritImageSet, foveaRegion );
    for( meritImageSetIterator.GoToBegin(); !meritImageSetIterator.IsAtEnd();
         ++meritImageSetIterator )
      {
      meritImage = meritImageSetIterator.Get();
      
      ImageRegionIteratorWithIndex<MeritImageType>
        meritImageIterator( meritImage, foveaRegion );
      for( meritImageIterator.GoToBegin(); !meritImageIterator.IsAtEnd();
           ++meritImageIterator )
        {
        trailMerit = meritImageIterator.Get();
        if(  ( trailMerit > bestTrailMerit )  &&  ( openPath || 
          (meritImageSetIterator.GetIndex()==meritImageIterator.GetIndex()) )  )
          {
          // this is the best trail we've seen so far
          bestTrailMerit = trailMerit;
          bestTrailStartFoveaIndex = meritImageSetIterator.GetIndex();
          bestTrailEndFoveaIndex   = meritImageIterator.GetIndex();
          }
        }
      }
    itkDebugMacro( << "GenerateData():  Best trail start fovea index = " << bestTrailStartFoveaIndex );
    
    
    // Convert the optimal trail into the optimal path
    itkDebugMacro( << "GenerateData():  Optimal trail -> path." );
    
    // Compute the new starting location for the optimal path
    OffsetType StartOfPathCorrection;
    for( unsigned int i=0; i<ImageDimension; i++ )
      StartOfPathCorrection[i] = bestTrailStartFoveaIndex[i] - m_FoveaRadius;
    outputPath->SetStart( inputPath->GetStart() + StartOfPathCorrection );
    
    // Get a single, consecutive list of the optimal steps to form the path.
    // We must work backwords...
    
    // Fill dirtyPath with the "dirty" optimal trail.
    bestTrailCurrentFoveaIndex = bestTrailEndFoveaIndex;
    OffsetType step;
    OffsetType equivalentStep;
    OffsetType zeroOffset;
    zeroOffset.Fill(0);
    // Extract the list of computed optimal steps:
    for( int i = pathLength - 1; i >= 0; i-- )
      {
      step = tstepImageSet[i]->GetPixel(bestTrailStartFoveaIndex)
                                  ->GetPixel(bestTrailCurrentFoveaIndex);
      bestTrailCurrentFoveaIndex += inputPath->Evaluate(i) - step;
      
      // There are 5 cases to consider, each requiring additional computation:
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      
      // Do we want to avoid trying to trim the path?
      if( m_DoNotTrimPath )
      {
        // Keep all steps.
        dirtyPath.push_front( step );
        continue;
      }
      
      // Are there other steps w/ which to combine?
      if( dirtyPath.empty() )
        {
        // Keep this step (for now).
        dirtyPath.push_front( step );
        continue;
        }
      
      // Do these two steps cancel?
      equivalentStep = step + dirtyPath.front();
      if( equivalentStep == zeroOffset )
        {
        // Discard both steps.
        dirtyPath.pop_front();
        continue;
        }
      
      // Can these two steps be combined?
      unsigned int dim; // current dimension we are evaluating
      for( dim=0;  ( dim<ImageDimension )  &&  ( vnl_math_abs(equivalentStep[dim]) <= 1 );  dim++ )
        {
        // Do nothing; just loop.
        }
      if( dim == ImageDimension )
        {
        // Combine these steps.
        dirtyPath.front() = equivalentStep;
        continue;
        }
      
      // Now we know that these two steps can not be combined.
      // Keep both steps (for now).
      dirtyPath.push_front( step );
      
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // End of the 4 cases.
      
      } // loop through optimal steps
    } // end of this code-block
  
  
  
  // ----------------------------------------------------------------------- //
  //   Clean-up the  optimal trail & convert to a chain code path
  // ----------------------------------------------------------------------- //
    {
    if( !m_DoNotTrimPath )
      {
      itkDebugMacro( << "GenerateData():  Cleaning up the path." );
      
      // Adjust the STL list "dirtyPath" to clean up the trail before saving it as our path.
      
      // For closed paths, begin by checking if the beginning and ending n steps of the path can be combined.
      if( !openPath )
        {
        
        }
      
      // If two consecutive steps only move one voxel away, then replace with their single equivalent step.
      // If two consecutive steps completely cancel, then drop both of them.
      for( OffsetListIteratorType it =  dirtyPath.begin();
                                  it != --dirtyPath.end();  it++ )
        {
           
        }
      
      // George's Procedure:
      // Convert 8-to-4 (simplifies redundant removal)
      // RemoveRedundant (call 3 times)
      // while( ! path->validate() ) path->Ligate()
      // make path clockwise
      
      // My idea:
      // Put path on stack for removing redundant steps, so as not
      // to require a fixed number of calls to RemoveRedundant().
      // Also, we could examing the merit of single pixel wide
      // intrusions/extrusions to determine whether or not they should
      // be preserved.  By default, all single pixel wide "bunny trails"
      // will be removed.
      // 
      // I need to look more closely at George's Ligate and Validate code
      // to know for sure how I want to augment and implement his ideas.
      // 
      // I don't know that there is any reason for this filter to make the
      // path clockwise.
      
      } // End of the path cleaning
      
    // Build the optimal path from the cleaned-up list of optimal steps
    pathLength=0;
    for( OffsetListIteratorType it=dirtyPath.begin(); it!=dirtyPath.end(); it++ )
    {
      outputPath->InsertStep( pathLength++, *it );
    }
    
    } // End of this code block.
  
  
  
  // ----------------------------------------------------------------------- //
  //   Delete lingering data structures
  // ----------------------------------------------------------------------- //
  
  itkDebugMacro( << "GenerateData():  Freeing manually allocated memory." );

  delete [] tstepImageSet;
  // all other pointers are smart pointers and should automatically self-delete
}



// ------ Private Functions -------

template <class TChainCodePath, class TImage>
void
SwathChainCodePathFilter<TChainCodePath, TImage>
::FindPstepNeighborhoodOffsets(OffsetListType &offsetList,const OffsetType step)
{
  // We need to find the tail locations of all psteps coresponding to step by
  // finding the offsets from the tail of step to the tails of all psteps.
  // 
  // The offset-dimensional-components for which step has a non-zero component
  // have one set of choices, and the components for which step has a zero
  // component have another set of choices.  Therefore, we find all possible
  // offsets if we restrict ourselves to dimensions crossed by step, and we also
  // find all possible offsets if we likewise restrict ourselves to dimensions
  // NOT crossed by step.  Then, we take every combination of offsets from the
  // first list with the offsets from the second list.
  
  OffsetListType traversedList, notTraversedList;
  OffsetType offset;
  OffsetListIteratorType it;
  
  
  // Build traversedList, which coresponds to non-zero components of step, and
  // nonTraversedList, which coresponds to zero components of step.
  
  // start each list with the zero-offset
  offset.Fill(0);
  traversedList.push_front(offset);
  notTraversedList.push_front(offset);
  
  // grow the list appropriate for each component of step
  for( unsigned int i=0; i<ImageDimension; i++ )
    {
    if( 0 != step[i] )
      {
      // add elements to traversedList
      for( it = traversedList.begin(); it != traversedList.end(); it++)
        {
        offset = *it;
        offset[i] = step[i];
        if( offset != step )  // This check is what requires 2 seperate lists
          {
          traversedList.insert(it, offset); // inserted BEFORE current position
          }
        }
      }
    else
      {
      // add elements to notTraversedList
      for( it = notTraversedList.begin(); it != notTraversedList.end(); it++)
        {
        offset = *it;
        offset[i] = -1;
        notTraversedList.insert(it, offset); // inserted BEFORE current position
        offset[i] = +1;
        notTraversedList.insert(it, offset); // inserted BEFORE current position
        }
      }
    }
  
  
  // "Merge" traversedList and notTraversedList, taking all combinations
  offsetList.clear();
  for( OffsetListIteratorType traversedIt = traversedList.begin();
       traversedIt != traversedList.end();   traversedIt++ )
    for( OffsetListIteratorType notTraversedIt = notTraversedList.begin();
         notTraversedIt != notTraversedList.end();   notTraversedIt++ )
      {
      offsetList.push_front( *traversedIt + *notTraversedIt );
      }
  
  // Done
}

} // end namespace itk

#endif
