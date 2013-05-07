//vector class


class vect{
	public:
		double x;
		double y;
		double z;
		vect(){
			x=0;
			y=0;
			z=0;
		}
		vect(double xin, double yin, double zin){
			x=xin;
			y=yin;
			z=zin;
		}
		bool operator==(const vect &rhs){
			if ((rhs.x==x)&&(rhs.y==y)&&(rhs.z==z)){
				return true;
			}
			return false;
		}
		vect & operator=(const vect &rhs){
			if (this==&rhs) return *this;
			this->x=rhs.x;
			this->y=rhs.y;
			this->z=rhs.z;
			return *this;
		}
		vect & operator+=(const vect &rhs){
			x=x+rhs.x;
			y=y+rhs.y;
			z=z+rhs.z;
			return *this;
		}
		vect & operator-=(const vect &rhs){
			x=x-rhs.x;
			y=y-rhs.y;
			z=z-rhs.z;
			return *this;
		}
		vect & operator*=(const double &rhs){
			x=x*rhs;
			y=y*rhs;
			z=z*rhs;
			return *this;
		}
		vect operator+(const vect &other){
			vect result=*this;
			result+=other;
			return result;
		}
		vect operator-(const vect &other){
			vect result=*this;
			result-=other;
			return result;
		}
		vect operator*(const double &other){
			vect result=*this;
			result*=other;
			return result;
		}
		bool operator!=(const vect& other){
			return !(*this==other);
		}
		double mag(){
			return sqrt(x*x+y*y+z*z);
		}
};


double dot(const vect& a, const vect& b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline vect cross(const vect& a, const vect& b){
     return vect(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

void print(vect a){
    cout<<"["<<a.x<<", "<<a.y<<", "<<a.z<<"]\n";
}