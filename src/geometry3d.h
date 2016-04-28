/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _geomety3d_H
#define _geomety3d_H
#include <cmath>
#include <vector>


class cVec;//forward declarations
class cPnt;//forward declarations
class cLine;//forward declarations
class cLineSeg;//forward declarations


class cVec{
protected:
public:

	double x;
	double y;
	double z;

	//cVec(const cVec& v){x=v.x;y=v.y;z=v.z;}
	cVec(){x=0.0; y=0.0; z=0.0;}
	cVec(double xo, double yo, double zo){x=xo; y=yo; z=zo;}	
	inline void set(double xo, double yo, double zo){x=xo; y=yo; z=zo;}
	
	int operator==(cVec a);
	int operator!=(cVec a);
	inline cVec& operator=(double a){x=a;y=a;z=a;return *this;}		
	inline cVec& operator=(cVec v){x=v.x;y=v.y;z=v.z;return *this;}		
	cVec& operator=(const std::vector<double>& v){
		x = v[0]; y = v[1]; z = v[2];
		return *this;
	}
	cVec(const std::vector<double>& v){
		x = v[0]; y = v[1]; z = v[2];
	}
	
	
	inline cVec& operator+=(cVec v){
		x += v.x; y += v.y; z += v.z;
		return *this;
	}
	inline cVec& operator-=(cVec v){
		x -= v.x; y -= v.y; z -= v.z;
		return *this;
	}		
	inline cVec& operator*=(cVec v){
		x *= v.x; y *= v.y; z *= v.z;
		return *this;
	}		
	inline cVec& operator/=(cVec v){
		x /= v.x; y /= v.y; z /= v.z;
		return *this;
	}		
	inline cVec& operator*=(double s){
		x*=s; y*=s; z*=s;
		return *this;
	}	
	inline cVec& operator/=(double s){
		x/=s; y/=s; z/=s;
		return *this;
	}
	
	

	inline cVec operator+(const cVec a){return cVec(x+a.x, y+a.y, z+a.z);}	
	inline cVec operator-(const cVec a){return cVec(x-a.x, y-a.y, z-a.z);}
	inline cVec operator*(cVec a){return cVec(x*a.x, y*a.y, z*a.z);}
	inline cVec operator/(cVec a){return cVec(x/a.x, y/a.y, z/a.z);}
	inline cVec operator-(){return cVec (-x, -y, -z);}	
	inline friend cVec operator+(const cVec& a, const cVec& b){return cVec(a.x+b.x, a.y+b.y, a.z+b.z);}
	inline friend cVec operator-(const cVec& a, const cVec& b){return cVec(a.x-b.x, a.y-b.y, a.z-b.z);}
	inline friend cVec operator*(cVec v, double s){return cVec(v.x*s, v.y*s, v.z*s);}
	inline friend cVec operator*(double s, cVec v){return cVec(v.x*s, v.y*s, v.z*s);}
	inline friend cVec operator/(cVec v, double s){return cVec(v.x/s, v.y/s, v.z/s);}
	inline friend cVec operator/(double s, cVec v){return cVec(v.x/s, v.y/s, v.z/s);}
	

	inline double length_squared(){
		return x*x + y*y + z*z;
	}
	inline double length(){
		return sqrt(x*x + y*y + z*z);
	}
	inline double length2n(double power)
	{
		return pow(x*x + y*y + z*z,0.5*power);
	}
	void scale(double scalefactor);
	cVec scaledcopy(double scalefactor);	
	inline void unitise(){	
		double len = sqrt(x*x + y*y + z*z);
		x/=len; y/=len; z/=len;
	}
	
	inline cVec unit()
	{
		cVec temp = *this;
		temp.unitise();
		return temp;
	}
	
	inline double dot(const cVec v) const
	{
		return v.x*x + v.y*y + v.z*z;
	}
	inline cVec cross(const cVec b) const 
	{		
		cVec c;//c = a x b (a is this cVec)
		c.x = y*b.z - b.y*z;
		c.y = z*b.x - b.z*x;
		c.z = x*b.y - b.x*y;
		return c;
	}		
	

	friend double distance(cVec a, cVec b);
	friend double distance2n(cVec a, cVec b, double power);
	friend cVec cross(cVec a, cVec b);
	friend double dot(cVec a, cVec b);
	friend cVec scale(cVec m, double scalefactor);
	friend cVec unit(cVec m){return m.unit();}

	cVec rotate(double angle, cVec axis);

};
class cPnt : public cVec{

public:
	
	cPnt(double xo, double yo, double zo){x=xo; y=yo; z=zo;}	
	cPnt(){x=0.0; y=0.0; z=0.0;}
				
	cPnt& operator=(const std::vector<double>& v){
		x = v[0]; y = v[1]; z = v[2];
		return *this;
	}
	cPnt(const std::vector<double>& v){
		x = v[0]; y = v[1]; z = v[2];
	}

	inline cPnt operator=(cVec a){return cPnt(a.x, a.y, a.z);}
	inline cPnt operator+(cVec a){return cPnt(x+a.x, y+a.y, z+a.z);}
	inline cPnt operator-(cVec a){return cPnt(x-a.x, y-a.y, z-a.z);}

	inline cPnt operator+(const cPnt a){return cPnt(x+a.x, y+a.y, z+a.z);}	
	inline cPnt operator-(const cPnt a){return cPnt(x-a.x, y-a.y, z-a.z);}
	inline cPnt operator*(cPnt a){return cPnt(x*a.x, y*a.y, z*a.z);}
	inline cPnt operator/(cPnt a){return cPnt(x/a.x, y/a.y, z/a.z);}
	inline cPnt operator-(){return cPnt (-x, -y, -z);}	
	inline friend cPnt operator+(const cPnt& a, const cPnt& b){return cPnt(a.x+b.x, a.y+b.y, a.z+b.z);}
	inline friend cPnt operator-(const cPnt& a, const cPnt& b){return cPnt(a.x-b.x, a.y-b.y, a.z-b.z);}
	inline friend cPnt operator*(cPnt v, double s){return cPnt(v.x*s, v.y*s, v.z*s);}
	inline friend cPnt operator*(double s, cPnt v){return cPnt(v.x*s, v.y*s, v.z*s);}
	inline friend cPnt operator/(cPnt v, double s){return cPnt(v.x/s, v.y/s, v.z/s);}
	inline friend cPnt operator/(double s, cPnt v){return cPnt(v.x/s, v.y/s, v.z/s);}

	//inline cPnt operator+(cPnt p){return cPnt(x+p.x, y+p.y, z+p.z);}
	//inline friend cPnt operator+(const cPnt& a, const cPnt& b){return cPnt(a.x+b.x, a.y+b.y, a.z+b.z);}
	
	double distance(cPnt p);
	int online(cLine m);
};
class cLine{

protected:

	cVec PLL;
	cPnt PNT;

public:

	cLine(cPnt p, cPnt q){PLL=unit(p-q); PNT=p;}
	cLine(cVec v, cPnt p){PLL=v.unit(); PNT=p;}
	cLine(){PLL=cVec(); PNT=cPnt();}

	inline cVec pll(){return PLL;}
	inline cPnt pnt(){return PNT;}

	inline void setpll(cVec v){PLL=v.unit();}
	inline void setpnt(cPnt p){PNT=p;}
	inline void set(cVec v, cPnt p){PLL=v.unit(); PNT=p;}

	int operator==(cLine m);
	int operator!=(cLine m);
	cLine operator=(cLine m);	

	friend double distance(cLine m, cPnt p);
	friend double distance(cLine m, cLine n);
	friend cPnt closestpoint(cLine m, cPnt p);	
	friend int closestpoints(cLine m, cLine n, cPnt &p, cPnt &q);
};

cPnt closestpoint(cLineSeg s, cPnt p);


class cLineSeg{

protected:

	cPnt P;
	cPnt Q;

public:

	cLineSeg(){P=cPnt(); Q=cPnt();}
	cLineSeg(cPnt p, cPnt q){P=p; Q=q;}

	inline cPnt p(){return P;}
	inline cPnt q(){return Q;}

	inline void setp(cPnt p){P=p;}
	inline void setq(cPnt q){Q=q;}
	inline void set(cPnt p, cPnt q){P=p; Q=q;}
};

inline cVec unitnormal(const cPnt& p1, const cPnt& p2, const cPnt& p3)
{		
	//p1 p2 p3 are counter clockwise
	//about the returned normal
	cVec a = p2-p1;
	cVec b = p3-p2;	
	cVec n = a.cross(b);	
	n.unitise();			
	return n;
}
#endif

