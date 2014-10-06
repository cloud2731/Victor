/*  This file is part of Victor.

    Victor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Victor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Victor.  If not, see <http://www.gnu.org/licenses/>.
*/
 

#ifndef __JensenShannon_H__
#define __JensenShannon_H__

#include <Profile.h>
#include <ScoringFunction.h>

namespace Biopool
{
/** @brief Calculate scores for profile to profile alignment using
*                  Jensen-Shannon method. 
 * 
* @Description  Some explanations can be found in:
*
*                  Guoli Wang, Roland L. Dunbrack jr.
*                  Scoring profile-to-profile sequence alignments.
*                  Institute for Cancer Research, Fox Chase Cancer Center,
*                  Philadelphia, Pennsylvania 19111, USA. March 16, 2004.
* @This 
 **/
class JensenShannon : public ScoringFunction
{

public:

// CONSTRUCTORS:

	/// Default constructor.
	JensenShannon(Profile *pro1, Profile *pro2);

	/// Copy constructor.
	JensenShannon(const JensenShannon &orig);

	/// Destructor.
	virtual ~JensenShannon();


// OPERATORS:

	/// Assignment operator.
	JensenShannon& operator = (const JensenShannon &orig);


// PREDICATES:

	/// Calculate scores to create matrix values.
	virtual double scoringSeq(int i, int j);


// MODIFIERS:

	/// Copy orig object to this object ("deep copy").
	virtual void copy(const JensenShannon &orig);

	/// Construct a new "deep copy" of this object.
	virtual JensenShannon* newCopy();


protected:


private:

// ATTRIBUTES:

	Profile *pro1;    ///< Target profile.
	Profile *pro2;    ///< Template profile.
	double p1[20];    ///< Target background frequencies.
	double p2[20];    ///< Template background frequencies.

};

} // namespace

#endif