/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkParametricBlindLeastSquaresDeconvolutionImageFilter_h
#define itkParametricBlindLeastSquaresDeconvolutionImageFilter_h

#include "itkIterativeDeconvolutionImageFilter.h"

#include "itkParametricImageSource.h"
#include "itkTernaryFunctorImageFilter.h"

namespace itk
{
/** \class ParametricBlindLeastSquaresDeconvolutionImageFilter
 *
 * \brief Least-squares blind deconvolution filter that also estimates
 * the parameters of a user-supplied parametric point-spread function.
 *
 * This filter takes a parametric kernel image source instead of a
 * static kernel image. During the deconvolution iterations, a new
 * estimate of the restored image is produced, along with a new
 * estimate of the kernel parameters. The parameters are available
 * through the kernel image sources GetParameters() method after the
 * filter has executed.
 *
 * Both the image estimate and the kernel parameter estimate are
 * produced through gradient descent on a sum-of-squared differences
 * objective function, making this method suitable for zero-mean
 * Gaussian white noise.
 *
 * This filter produces output in two forms: a deconvolved image and the
 * parameters of the input kernel image source.
 *
 * \warning The method SetKernelImage(), inherited from the superclass
 * of this filter, is publicly available. However, this algorithm does
 * not use the static kernel image set through this method. Instead,
 * it uses the output of the parametric kernel source you specify.
 *
 * \author Cory Quammen, The University of North Carolina at Chapel Hill
 *
 * \ingroup ITKDeconvolution
 */
template< typename TInputImage, typename TKernelSource, typename TOutputImage=TInputImage >
class ITK_TEMPLATE_EXPORT ParametricBlindLeastSquaresDeconvolutionImageFilter :
    public IterativeDeconvolutionImageFilter< TInputImage, typename TKernelSource::OutputImageType, TOutputImage >
{
public:
  ITK_DISALLOW_COPY_AND_ASSIGN(ParametricBlindLeastSquaresDeconvolutionImageFilter);

  /** Standard type alias. */
  using Self = ParametricBlindLeastSquaresDeconvolutionImageFilter;
  using Superclass = IterativeDeconvolutionImageFilter< TInputImage,
                                             typename TKernelSource::OutputImageType,
                                             TOutputImage >;
  using Pointer = SmartPointer< Self >;
  using ConstPointer = SmartPointer< const Self >;

  /** Other useful type alias. */
  using InputImageType = TInputImage;
  using OutputImageType = TOutputImage;

  /** Internal types used by the FFT filters. */
  using InternalImageType = typename Superclass::InternalImageType;
  using InternalImagePointerType = typename Superclass::InternalImagePointerType;
  using InternalComplexType = typename Superclass::InternalComplexType;
  using InternalComplexImageType = typename Superclass::InternalComplexImageType;
  using InternalComplexImagePointerType = typename Superclass::InternalComplexImagePointerType;

  /** Type for the parametric kernel source. */
  using KernelSourceType = TKernelSource;
  using KernelSourcePointer = typename KernelSourceType::Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Runtime information support. */
  itkTypeMacro(ParametricBlindLeastSquaresDeconvolutionImageFilter,
               IterativeDeconvolutionImageFilter);

  /** Set/get the parametric kernel source. */
  void SetKernelSource(KernelSourceType * kernelSource);
  itkGetModifiableObjectMacro(KernelSource, KernelSourceType);

  /** Set/get the scale factor (also known as learning rate) for the
   * image intensity gradient descent. */
  itkSetMacro(Alpha, double);
  itkGetMacro(Alpha, double);

  /** Set/get the scale factor (also known as learning rate) for the
   * parameter gradient descent. */
  itkSetMacro(Beta, double);
  itkGetMacro(Beta, double);

protected:
  ParametricBlindLeastSquaresDeconvolutionImageFilter();
  ~ParametricBlindLeastSquaresDeconvolutionImageFilter() override = default;

  void Initialize(ProgressAccumulator * progress,
                          float progressWeight,
                          float iterationProgressWeight) override;

  void Iteration(ProgressAccumulator * progress,
                         float iterationProgressWeight) override;

  void Finish(ProgressAccumulator *progress, float progressWeight) override;

  void PrintSelf(std::ostream & os, Indent indent) const override;

private:
  template< typename TPixel >
    class ITK_TEMPLATE_EXPORT ParametricBlindLeastSquaresDeconvolutionDifference
  {
  public:
    ParametricBlindLeastSquaresDeconvolutionDifference() = default;
    ~ParametricBlindLeastSquaresDeconvolutionDifference() = default;

    bool operator!=(const ParametricBlindLeastSquaresDeconvolutionDifference &) const
    {
      return false;
    }

    bool operator==(const ParametricBlindLeastSquaresDeconvolutionDifference & other) const
    {
      return !( *this != other );
    }

    inline TPixel operator()(const TPixel & estimateFT,
                             const TPixel & kernelEstimateFT,
                             const TPixel & inputFT)
    {
      return estimateFT * kernelEstimateFT - inputFT;
    }
  };

  template< typename TPixel >
    class ITK_TEMPLATE_EXPORT ParametricBlindLeastSquaresDeconvolutionImageUpdate
  {
  public:
    ParametricBlindLeastSquaresDeconvolutionImageUpdate()  {}
    ~ParametricBlindLeastSquaresDeconvolutionImageUpdate() = default;

    bool operator!=(const ParametricBlindLeastSquaresDeconvolutionImageUpdate &) const
    {
      return false;
    }

    bool operator==(const ParametricBlindLeastSquaresDeconvolutionImageUpdate & other) const
    {
      return !( *this != other );
    }

    inline TPixel operator()(const TPixel & estimateFT,
                             const TPixel & differenceFT,
                             const TPixel & kernelFT)
    {
      // Because of the linearity of the Fourier transform, we can
      // perform the update step in the Fourier domain
      return estimateFT - m_Alpha * ( differenceFT * std::conj( kernelFT ) );
    }

    void SetAlpha(double alpha)
    {
      m_Alpha = alpha;
    }
    double GetAlpha() const
    {
      return m_Alpha;
    }

  private:
    double m_Alpha{0.01};
  };

  KernelSourcePointer             m_KernelSource;

  /** Step sizes for the gradient descent of the image and the
   * kernel parameters. These are very different spaces, so they
   * deserve different step size parameters. */
  double m_Alpha;
  double m_Beta;

  /** Temporary images. */
  InternalComplexImagePointerType m_TransformedInput;
  InternalComplexImagePointerType m_TransformedCurrentEstimate;

  /** These are the internal filters that perform the updating of the
   * image estimate. */
  using DifferenceFunctorType = ParametricBlindLeastSquaresDeconvolutionDifference<InternalComplexType>;
  using DifferenceFilterType = TernaryFunctorImageFilter< InternalComplexImageType,
    InternalComplexImageType,
    InternalComplexImageType,
    InternalComplexImageType,
    DifferenceFunctorType >;
  typename DifferenceFilterType::Pointer m_DifferenceFilter;

  using ImageUpdateFunctorType = ParametricBlindLeastSquaresDeconvolutionImageUpdate<InternalComplexType>;
  using ImageUpdateFilterType = TernaryFunctorImageFilter< InternalComplexImageType,
    InternalComplexImageType,
    InternalComplexImageType,
    InternalComplexImageType,
    ImageUpdateFunctorType >;
  typename ImageUpdateFilterType::Pointer m_ImageUpdateFilter;

};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkParametricBlindLeastSquaresDeconvolutionImageFilter.hxx"
#endif


#endif
