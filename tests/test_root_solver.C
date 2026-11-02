#include "../root/root_solver.hh"
#include <iostream>
namespace {
void close_matrix(const TMatrixD& a,const TMatrixD& b){
    for(int i=0;i<a.GetNrows();i++)for(int j=0;j<a.GetNcols();j++)
        if(std::abs(a(i,j)-b(i,j))>1e-8*(1+std::abs(b(i,j))))throw std::runtime_error("ROOT solver invariant failed");
}
}
void test_root_solver(){
    using namespace root_deconv;const int n=6,d=7,dim=14;
    CenteredScore constant,nearConstant,withZeros;
    std::vector<double> nearScores;
    for(int i=0;i<50;i++){
        constant.add(1./50);
        double x=(1.+(i%2?1e-7:-1e-7))/50;nearConstant.add(x);nearScores.push_back(x);
        if(i<10)withZeros.add(1./50);
    }
    if(constant.variance_of_sum(50)!=0)throw std::runtime_error("Constant-score variance is not exactly zero");
    long double mean=0,oracle=0;for(double x:nearScores)mean+=x;mean/=50;
    for(double x:nearScores)oracle+=(x-mean)*(x-mean);oracle*=50.L/49;
    if(std::abs(nearConstant.variance_of_sum(50)/oracle-1)>1e-9)throw std::runtime_error("Small real variance was lost");
    if(std::abs(withZeros.variance_of_sum(50)-.0032*50/49)>1e-15)throw std::runtime_error("Implicit zero histories were lost");
    CenteredScore sparseX,sparseY;CenteredPair sparsePair;
    std::vector<std::pair<double,double>> sparse;
    long double mx=0,my=0,expectedCross=0;
    for(int i=0;i<257;i++){
        double x=i%7?100000.+(i%3)*.01:0,y=i%5?-10*x+(i%13)*.005:0;
        sparse.emplace_back(x,y);sparseX.add(x);sparseY.add(y);sparsePair.add(x,y);mx+=x;my+=y;
    }
    mx/=sparse.size();my/=sparse.size();
    for(auto xy:sparse)expectedCross+=(xy.first-mx)*(xy.second-my);expectedCross*=257.L/256;
    if(std::abs(sparsePair.covariance_of_sum(sparseX,sparseY,257)/expectedCross-1)>1e-12)
        throw std::runtime_error("Sparse pair covariance disagrees with independent centered sum");
    TMatrixD y(5,dim),v(5*dim,5*dim),b(n,n);
    for(int c=0;c<5;c++){
        double ref=20*(c+1);y(c,n)=1e7;y(c,d+n)=1e7*ref;
        for(int i=0;i<n;i++){
            y(c,i)=1e6*(c==i?5:1)*(i+1);y(c,d+i)=y(c,i)*ref*(1+.01*(i-c));b(i,i)+=y(c,i);
        }
    }
    for(int i=0;i<5*dim;i++)v(i,i)=1e7;
    auto fit=solve(y,v,b);
    TMatrixD independent(5,5);
    for(int i=0;i<n;i++){
        double rate=0;for(int c=0;c<5;c++)rate+=y(c,i);
        if(std::abs(fit.C(i,i)*(.8*.8*344*86400*rate)/1e18-1)>1e-12)throw std::runtime_error("Independent covariance limit failed");
        for(int c=0;c<5;c++)for(int e=0;e<5;e++)independent(c,e)+=fit.H(i,c)*fit.H(i,e)/fit.C(i,i);
    }
    close_matrix(fit.M,independent);
    // Change signal units separately in every detector, leaving reference counts.
    TMatrixD gy=y,gv=v,gb=b;double gain[]={.01,2.,.1,10.,3.,100.};
    for(int c=0;c<5;c++)for(int i=0;i<n;i++){
        gy(c,i)*=gain[i];gy(c,d+i)*=gain[i];
        gv(c*dim+i,c*dim+i)*=gain[i]*gain[i];gv(c*dim+d+i,c*dim+d+i)*=gain[i]*gain[i];
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)gb(i,j)*=gain[i]*gain[j];
    auto scaled=solve(gy,gv,gb);
    close_matrix(fit.C,scaled.C);close_matrix(fit.H,scaled.H);close_matrix(fit.M,scaled.M);
    close_matrix(fit.experimental_covariance,scaled.experimental_covariance);close_matrix(fit.mc_covariance,scaled.mc_covariance);
    // Deliberately make all component columns identical, so they cannot be fitted.
    TMatrixD degenerate=y;
    for(int c=0;c<5;c++)for(int i=0;i<n;i++){
        degenerate(c,i)=y(0,i);degenerate(c,d+i)=y(0,i)*20*(c+1);
    }
    bool rejected=false;try{solve(degenerate,v,b);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Rank deficiency was not rejected");
    TMatrixD invalid=v;invalid(0,0)=std::numeric_limits<double>::quiet_NaN();
    rejected=false;try{solve(y,invalid,b);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Nonfinite MC input was not rejected");
    invalid.ResizeTo(b);invalid=b;invalid(0,1)=1;
    rejected=false;try{solve(y,v,invalid);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Asymmetric B was not rejected");
    invalid.ResizeTo(v);invalid=v;invalid(0,1)=invalid(1,0)=2e7;
    rejected=false;try{solve(y,invalid,b);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Indefinite MC covariance was not rejected");
    invalid=v;invalid(0,0)=-1e-20;
    rejected=false;try{solve(y,invalid,b);}catch(const std::runtime_error&){rejected=true;}
    if(!rejected)throw std::runtime_error("Negative MC covariance was silently clamped");
    // A false global rate normalization has a known, nonzero effect on the fit.
    TMatrixD twiceY=y,twiceB=b,fourV=v;twiceY*=2;twiceB*=2;fourV*=4;
    auto twice=solve(twiceY,fourV,twiceB);
    close_matrix(fit.F,twice.F);
    TMatrixD half=fit.experimental_covariance;half*=.5;
    close_matrix(half,twice.experimental_covariance);close_matrix(fit.mc_covariance,twice.mc_covariance);
    // Closed-form oracle: five pure detectors, H=I, hence Cov(A)=C exactly.
    TMatrixD pureY(5,12),pureV(60,60),pureB(5,5),knownMC(5,5),knownM(5,5);
    for(int c=0;c<5;c++){
        pureY(c,c)=1e6;pureY(c,6+c)=1e6*10*(c+1);
        pureY(c,5)=100;pureY(c,11)=100*10*(c+1);
        pureB(c,c)=(c+1)*(c+1);knownM(c,c)=1/pureB(c,c);
        pureV(c*12+11,c*12+11)=100*(c+1)*(c+1);
        knownMC(c,c)=.01*(c+1)*(c+1);
    }
    auto pure=solve(pureY,pureV,pureB,1e6/86400,1);
    close_matrix(pure.experimental_covariance,pureB);close_matrix(pure.M,knownM);
    close_matrix(pure.mc_covariance,knownMC);
    // Additional background: six pure channels must produce a six-component fit.
    TMatrixD sixY(6,14),sixV(84,84),sixB(6,6);
    for(int c=0;c<6;c++){
        sixY(c,c)=1e6;sixY(c,7+c)=1e6*10*(c+1);
        sixY(c,6)=100;sixY(c,13)=100*10*(c+1);sixB(c,c)=(c+1)*(c+1);
    }
    auto six=solve(sixY,sixV,sixB,1e6/86400,1);
    if(six.M.GetNrows()!=6)throw std::runtime_error("Additional background dimension lost");
    close_matrix(six.experimental_covariance,sixB);
    std::cout<<"PASS: six-component background extension matches closed-form covariance\n";
    std::cout<<"PASS: ROOT independent limit, gain invariance, singular/invalid input rejection and false-rate scaling\n";
    std::cout<<"PASS: exact statistical sigmas 1,2,3,4,5 ppb; exact MC sigmas 0.1,0.2,0.3,0.4,0.5 ppb\n";
    std::cout<<"PASS: stable constant/small variances and implicit zero histories; invalid covariance still rejected\n";
}
