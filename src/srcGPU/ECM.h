class SceECM {
	SceNodes* nodes;
        void applyECMForce_M ; 
}


struct AddECMForce
{
    const float a;

    AddECMForce(float _a) : a(_a) {}

    __host__ __device__
        double operator()(const double & Y) const {

 if (Y<a){
            return a ;
        }
 else {     return Y ; 

}
};



