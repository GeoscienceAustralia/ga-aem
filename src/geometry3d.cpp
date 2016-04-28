/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cfloat>
#include <climits>

#include "general_constants.h"
#include "geometry3d.h"


//cVec///////////////////////////////////////
int cVec::operator==(cVec a)
{
	if(fabs(x - a.x)<DBL_EPSILON && fabs(y - a.y)<DBL_EPSILON && fabs(z - a.z)<DBL_EPSILON)return 1;
	else return 0;
}
int cVec::operator!=(cVec a)
{
	if (*this == a)return 0;
	else return 1;
}
void cVec::scale(double scalefactor)
{
	x=x*scalefactor;
	y=y*scalefactor;
	z=z*scalefactor;
}
cVec cVec::scaledcopy(double scalefactor)
{
	cVec temp(x*scalefactor,y*scalefactor,z*scalefactor);
	return temp;
}
cVec cVec::rotate(double angle, cVec axis)
{
	cVec a=axis.unit();
	double cosangle = cos(angle*D2R);
	double sinangle = sin(angle*D2R);
	
	double xx = a.x*a.x;   double yy = a.y*a.y;   double zz = a.z*a.z;
	double xy = a.x*a.y;   double xz = a.x*a.z;   double yz = a.y*a.z;

	double mat[9];

	double oneminuscosangle = (1.0-cosangle);

	mat[0]=xx+cosangle*(1-xx);
	mat[4]=yy+cosangle*(1-yy);
	mat[8]=zz+cosangle*(1-zz);

	mat[1]=xy*oneminuscosangle - a.z*sinangle;
	mat[3]=xy*oneminuscosangle + a.z*sinangle;

	mat[2]=xz*oneminuscosangle + a.y*sinangle;
	mat[6]=xz*oneminuscosangle - a.y*sinangle;

	mat[5]=yz*oneminuscosangle - a.x*sinangle;
	mat[7]=yz*oneminuscosangle + a.x*sinangle;

	return cVec(
		mat[0]*x + mat[1]*y + mat[2]*z,
		mat[3]*x + mat[4]*y + mat[5]*z,
		mat[6]*x + mat[7]*y + mat[8]*z
	);
}
//cPnt///////////////////////////////////////
int cPnt::online(cLine m)
{
	if(*this==m.pnt())return 1;
	cVec a=*this-m.pnt();
	a.unitise();
	if(m.pll()==a || m.pll()==-a)return 1;
	else return 0;
}
double cPnt::distance(cPnt p)
{
	return pow(pow(p.x-x,2.0) + pow(p.y-y,2.0) + pow(p.z-z,2.0) , 0.5);
}
//cLine///////////////////////////////////////
int cLine::operator==(cLine m)
{
	if(PLL==m.PLL){
		cVec v = PNT - m.PNT;
		v.unitise();
		if(PLL==v || PLL==-v)return 1;
	}
	return 0;
}
int cLine::operator!=(cLine m)
{
	if(*this==m) return 0;
	else return 1;
}
inline cLine cLine::operator=(cLine m)
{
	PLL=m.PLL;
	PNT=m.PNT;
	return *this;
}
//Friend Functions///////////////////////////////////////
inline double distance(cVec a, cVec b)
{
	return pow(pow(a.x-b.x,2.0) + pow(a.y-b.y,2.0) + pow(a.z-b.z,2.0) , 0.5);
}
inline double distance2n(cVec a, cVec b, double power)
{
	return pow(pow(a.x-b.x,2.0) + pow(a.y-b.y,2.0) + pow(a.z-b.z,2.0) , 0.5*power);
}
cVec cross(cVec a, cVec b)
{
	cVec c;
	c.x = a.y*b.z - b.y*a.z;
	c.y = a.z*b.x - b.z*a.x;
	c.z = a.x*b.y - b.x*a.y;
	return c;
}
double dot(cVec a, cVec b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
cVec scale(cVec m, double scalefactor)
{
	cVec temp;
	temp.x=m.x*scalefactor;
	temp.y=m.y*scalefactor;
	temp.z=m.z*scalefactor;
	return temp;
}
double distance(cLine m, cPnt p)
{
	cVec v=cross(m.PLL, (p-m.PNT));
	return v.length();
}
double distance(cLine m, cLine n)
{
	cVec w=m.PNT - n.PNT;
	if(w.length()<DBL_EPSILON)return 0;
	cVec v=cross(n.PLL, m.PLL);
	if(v.length()<DBL_EPSILON)return distance(n, m.PNT);
	return fabs(dot(w,v))/v.length();
}
cPnt closestpoint(cLine m, cPnt p)
{
	double s = dot(m.PLL, (p-m.PNT));
	return m.PNT + scale(m.PLL,s);
}
cPnt closestpoint(cLineSeg s, cPnt p)
{
	cLine m(s.p(),s.q());
	cPnt c = closestpoint(m,p);

	double dx = fabs(s.p().x - s.q().x);
	double dy = fabs(s.p().y - s.q().y);

	if(dx>dy){
		if(c.x<s.p().x && c.x<s.q().x){}
		else if(c.x>s.p().x && c.x>s.q().x){}
		else return c;
	}
	else if(dx<dy){
		if(c.y<s.p().y && c.y<s.q().y){}
		else if(c.y>s.p().y && c.y>s.q().y){}
		else return c;
	}

	double dp=distance(p,s.p());
	double dq=distance(p,s.q());

	if(dp < dq)return s.p();
	else return s.q();
}
int closestpoints(cLine m, cLine n, cPnt &p, cPnt &q)
{
	double AD=dot(m.PNT,n.PLL);
	double AB=dot(m.PNT,m.PLL);
	double BC=dot(m.PLL,n.PNT);
	double BD=dot(m.PLL,n.PLL);
	double CD=dot(n.PNT,n.PLL);

	if((BD*BD-1.0)<DBL_EPSILON)return 0;

	double t = (AD - AB*BD + BC*BD - CD)/(1.0-BD*BD);
	q = n.PNT + scale(n.PLL,t);

	double s = -1.0 * dot(m.PLL, (m.PNT - q));
	p = m.PNT + scale(m.PLL,s);

	return 1;
}


