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
#ifndef itkScaleSkewVersor3DTransform_hxx
#define itkScaleSkewVersor3DTransform_hxx

#include "itkScaleSkewVersor3DTransform.h"
#include "itkMath.h"

namespace itk
{
// Constructor with default arguments
template<typename TParametersValueType>
ScaleSkewVersor3DTransform<TParametersValueType>
::ScaleSkewVersor3DTransform() :
  Superclass(ParametersDimension)
{
  m_Scale.Fill(NumericTraits<TParametersValueType>::OneValue());
  m_Skew.Fill(NumericTraits<TParametersValueType>::ZeroValue());
}

// Constructor with arguments
template<typename TParametersValueType>
ScaleSkewVersor3DTransform<TParametersValueType>
::ScaleSkewVersor3DTransform(unsigned int parametersDimension) :
  Superclass(parametersDimension)
{
  m_Scale.Fill(1.0);
  m_Skew.Fill(0.0);
}

// Constructor with arguments
template<typename TParametersValueType>
ScaleSkewVersor3DTransform<TParametersValueType>::ScaleSkewVersor3DTransform(const MatrixType & matrix,
                                                                    const OutputVectorType & offset) :
  Superclass(matrix, offset)
{
  this->ComputeMatrixParameters();
}

// Directly set the matrix
template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetMatrix(const MatrixType & matrix)
{
  // Any matrix should work - bypass orthogonality testing
  using Baseclass = MatrixOffsetTransformBase<TParametersValueType, 3, 3>;
  this->Baseclass::SetMatrix(matrix);
}

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetMatrix(const MatrixType & matrix, const TParametersValueType itkNotUsed( tolerance ))
{
  // Any matrix should work - bypass orthogonality testing
  using Baseclass = MatrixOffsetTransformBase<TParametersValueType, 3, 3>;
  this->Baseclass::SetMatrix(matrix);
}

// Set Parameters
template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetParameters(const ParametersType & parameters)
{
  itkDebugMacro(<< "Setting parameters " << parameters);

  // Save parameters. Needed for proper operation of TransformUpdateParameters.
  if( &parameters != &(this->m_Parameters) )
    {
    this->m_Parameters = parameters;
    }

  // Transfer the versor part

  AxisType axis;

  double norm = parameters[0] * parameters[0];
  axis[0] = parameters[0];
  norm += parameters[1] * parameters[1];
  axis[1] = parameters[1];
  norm += parameters[2] * parameters[2];
  axis[2] = parameters[2];
  if( norm > 0 )
    {
    norm = std::sqrt(norm);
    }

  double epsilon = 1e-10;
  if( norm >= 1.0 - epsilon )
    {
    axis = axis / ( norm + epsilon * norm );
    }
  VersorType newVersor;
  newVersor.Set(axis);
  this->SetVarVersor(newVersor);

  itkDebugMacro(<< "Versor is now " << newVersor);

  // Matrix must be defined before translation so that offset can be computed
  // from translation
  m_Scale[0] = parameters[6];
  m_Scale[1] = parameters[7];
  m_Scale[2] = parameters[8];

  m_Skew[0] = parameters[9];
  m_Skew[1] = parameters[10];
  m_Skew[2] = parameters[11];
  m_Skew[3] = parameters[12];
  m_Skew[4] = parameters[13];
  m_Skew[5] = parameters[14];

  // Transfer the translation part
  TranslationType newTranslation;
  newTranslation[0] = parameters[3];
  newTranslation[1] = parameters[4];
  newTranslation[2] = parameters[5];

  this->SetVarTranslation(newTranslation);
  this->ComputeMatrix();
  this->ComputeOffset();

  // Modified is always called since we just have a pointer to the
  // parameters and cannot know if the parameters have changed.
  this->Modified();

  itkDebugMacro(<< "After setting parameters ");
}

//
// Get Parameters
//
// Parameters are ordered as:
//
// p[0:2] = right part of the versor (axis times std::sin(t/2))
// p[3:5] = translation components
// p[6:8] = Scale
// p[9:14] = Skew {xy, xz, yx, yz, zx, zy}
//

template<typename TParametersValueType>
const typename ScaleSkewVersor3DTransform<TParametersValueType>::ParametersType
& ScaleSkewVersor3DTransform<TParametersValueType>
::GetParameters() const
  {
  itkDebugMacro(<< "Getting parameters ");

  this->m_Parameters[0] = this->GetVersor().GetX();
  this->m_Parameters[1] = this->GetVersor().GetY();
  this->m_Parameters[2] = this->GetVersor().GetZ();

  this->m_Parameters[3] = this->GetTranslation()[0];
  this->m_Parameters[4] = this->GetTranslation()[1];
  this->m_Parameters[5] = this->GetTranslation()[2];

  this->m_Parameters[6] = this->GetScale()[0];
  this->m_Parameters[7] = this->GetScale()[1];
  this->m_Parameters[8] = this->GetScale()[2];

  this->m_Parameters[9] = this->GetSkew()[0];
  this->m_Parameters[10] = this->GetSkew()[1];
  this->m_Parameters[11] = this->GetSkew()[2];
  this->m_Parameters[12] = this->GetSkew()[3];
  this->m_Parameters[13] = this->GetSkew()[4];
  this->m_Parameters[14] = this->GetSkew()[5];

  itkDebugMacro(<< "After getting parameters " << this->m_Parameters);

  return this->m_Parameters;
  }

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetIdentity()
{
  m_Scale.Fill(NumericTraits<ScaleVectorValueType>::OneValue());
  m_Skew.Fill(NumericTraits<SkewVectorValueType>::ZeroValue());
  Superclass::SetIdentity();
}

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetScale(const ScaleVectorType & scale)
{
  m_Scale = scale;
  this->ComputeMatrix();
}

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::SetSkew(const SkewVectorType & skew)
{
  m_Skew = skew;
  this->ComputeMatrix();
}

// Compute the matrix
template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::ComputeMatrix()
{
  this->Superclass::ComputeMatrix();

  MatrixType newMatrix = this->GetMatrix();

  newMatrix[0][0] += m_Scale[0] - 1.0;
  newMatrix[1][1] += m_Scale[1] - 1.0;
  newMatrix[2][2] += m_Scale[2] - 1.0;
  newMatrix[0][1] += m_Skew[0];
  newMatrix[0][2] += m_Skew[1];
  newMatrix[1][0] += m_Skew[2];
  newMatrix[1][2] += m_Skew[3];
  newMatrix[2][0] += m_Skew[4];
  newMatrix[2][1] += m_Skew[5];

  this->SetVarMatrix(newMatrix);
}

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::ComputeMatrixParameters()
{
  itkExceptionMacro(<< "Setting the matrix of a ScaleSkewVersor3D transform is not supported at this time.");
}

// Print self
template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Scale:       " << m_Scale        << std::endl;
  os << indent << "Skew:        " << m_Skew         << std::endl;
}

template<typename TParametersValueType>
void
ScaleSkewVersor3DTransform<TParametersValueType>
::ComputeJacobianWithRespectToParameters(const InputPointType & p, JacobianType & jacobian) const
{
  using ValueType = typename VersorType::ValueType;

  // compute derivatives with respect to rotation
  const ValueType vx = this->GetVersor().GetX();
  const ValueType vy = this->GetVersor().GetY();
  const ValueType vz = this->GetVersor().GetZ();
  const ValueType vw = this->GetVersor().GetW();

  jacobian.SetSize( 3, this->GetNumberOfLocalParameters() );
  jacobian.Fill(0.0);

  const double px = p[0] - this->GetCenter()[0];
  const double py = p[1] - this->GetCenter()[1];
  const double pz = p[2] - this->GetCenter()[2];

  const double vxx = vx * vx;
  const double vyy = vy * vy;
  const double vzz = vz * vz;
  const double vww = vw * vw;

  const double vxy = vx * vy;
  const double vxz = vx * vz;
  const double vxw = vx * vw;

  const double vyz = vy * vz;
  const double vyw = vy * vw;

  const double vzw = vz * vw;

  // compute Jacobian with respect to quaternion parameters
  jacobian[0][0] = 2.0 * ( ( vyw + vxz ) * py + ( vzw - vxy ) * pz )
    / vw;
  jacobian[1][0] = 2.0 * ( ( vyw - vxz ) * px   - 2 * vxw   * py + ( vxx - vww ) * pz )
    / vw;
  jacobian[2][0] = 2.0 * ( ( vzw + vxy ) * px + ( vww - vxx ) * py   - 2 * vxw   * pz )
    / vw;

  jacobian[0][1] = 2.0 * ( -2 * vyw  * px + ( vxw + vyz ) * py + ( vww - vyy ) * pz )
    / vw;
  jacobian[1][1] = 2.0 * ( ( vxw - vyz ) * px                + ( vzw + vxy ) * pz )
    / vw;
  jacobian[2][1] = 2.0 * ( ( vyy - vww ) * px + ( vzw - vxy ) * py   - 2 * vyw   * pz )
    / vw;

  jacobian[0][2] = 2.0 * ( -2 * vzw  * px + ( vzz - vww ) * py + ( vxw - vyz ) * pz )
    / vw;
  jacobian[1][2] = 2.0 * ( ( vww - vzz ) * px   - 2 * vzw   * py + ( vyw + vxz ) * pz )
    / vw;
  jacobian[2][2] = 2.0 * ( ( vxw + vyz ) * px + ( vyw - vxz ) * py )
    / vw;

  jacobian[0][3] = 1.0;
  jacobian[1][4] = 1.0;
  jacobian[2][5] = 1.0;

  jacobian[0][6] = px;
  jacobian[1][7] = py;
  jacobian[2][8] = pz;

  jacobian[0][9]  = py;
  jacobian[0][10] = pz;
  jacobian[1][11] = px;
  jacobian[1][12] = pz;
  jacobian[2][13] = px;
  jacobian[2][14] = py;
}

} // namespace

#endif
