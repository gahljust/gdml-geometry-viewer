#pragma once
// ROOT-only equivalent of the launcher's validated GLS solver.
// Inputs: Y=[rate, rate*asymmetry] per process, joint finite-MC covariance V,
// and the physical second moment B (one rate weight, not rate squared).
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompChol.h>
#include <TDecompSVD.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
namespace root_deconv {
struct Fit {
    TMatrixD rates,F,F_mc_sigma,H,B,C,correlation,M;
    TMatrixD experimental_covariance,mc_covariance,combined_covariance;
    TVectorD reference_asymmetry_ppb;
    double scaled_condition=0;
};
inline TMatrixD transpose(const TMatrixD& m){return TMatrixD(TMatrixD::kTransposed,m);}
inline Fit solve(const TMatrixD& y,const TMatrixD& v,const TMatrixD& b,
                 double days=344,double polarization=.8){
    const int n=b.GetNrows(),d=n+1,dim=2*d,components=y.GetNrows();
    if(components<1||n<components||y.GetNcols()!=dim||v.GetNrows()!=components*dim||v.GetNcols()!=components*dim||b.GetNcols()!=n)
        throw std::runtime_error("Incompatible moment dimensions");
    if(!(std::isfinite(days)&&days>0&&std::isfinite(polarization)&&polarization>0&&polarization<=1))throw std::runtime_error("Invalid exposure or polarization");
    for(int c=0;c<components;c++)for(int i=0;i<dim;i++){
        if(!std::isfinite(y(c,i)))throw std::runtime_error("Nonfinite rate/asymmetry moment");
        if(i<d&&y(c,i)<0)throw std::runtime_error("Negative detector rate");
        if(i<d&&y(c,i)==0&&y(c,d+i)!=0)throw std::runtime_error("Asymmetry sum without detector rate");
    }
    auto validate_symmetric=[](const TMatrixD& m,const std::string& name,bool nonnegative){
        for(int i=0;i<m.GetNrows();i++)for(int j=0;j<m.GetNcols();j++){
            if(!std::isfinite(m(i,j)))throw std::runtime_error("Nonfinite "+name);
            if((nonnegative||i==j)&&m(i,j)<0)throw std::runtime_error("Negative "+name);
            if(std::abs(m(i,j)-m(j,i))>1e-10*std::max({1.,std::abs(m(i,j)),std::abs(m(j,i))}))
                throw std::runtime_error("Asymmetric "+name);
        }
    };
    validate_symmetric(b,"physical second moment",true);
    validate_symmetric(v,"MC covariance",false);
    // Allow singular/zero MC covariance, but reject negative-variance directions.
    // Normalize heterogeneous rate/asymmetry units; jitter is for validation only.
    std::vector<int> active;
    for(int i=0;i<v.GetNrows();i++){
        if(v(i,i)>0)active.push_back(i);
        else for(int j=0;j<v.GetNcols();j++)if(v(i,j)!=0)
            throw std::runtime_error("MC covariance with zero variance and nonzero covariance");
    }
    if(!active.empty()){
        int size=active.size();TMatrixD normalized(size,size);
        double tolerance=128*std::numeric_limits<double>::epsilon()*size;
        for(int i=0;i<size;i++)for(int j=0;j<size;j++){
            int a=active[i],b=active[j];double scale=std::sqrt(v(a,a))*std::sqrt(v(b,b));
            double x=v(a,b)/scale,z=v(b,a)/scale;
            if(!std::isfinite(x)||std::abs(x-z)>1e-10*std::max({1.,std::abs(x),std::abs(z)}))
                throw std::runtime_error("Asymmetric/nonfinite normalized MC covariance");
            normalized(i,j)=.5*(x+z)+(i==j?tolerance:0);
        }
        TDecompChol check(normalized);
        if(!check.Decompose())throw std::runtime_error("MC covariance is not positive semidefinite");
    }
    Fit f;f.rates.ResizeTo(n,components);f.F.ResizeTo(n,components);f.H.ResizeTo(n,components);f.B.ResizeTo(b);f.B=b;
    f.C.ResizeTo(n,n);f.correlation.ResizeTo(n,n);f.reference_asymmetry_ppb.ResizeTo(components);
    TVectorD r(n),sd(n),observed(n);
    for(int c=0;c<components;c++){
        if(!(y(c,n)>0))throw std::runtime_error("Missing Ring 5 open reference for component "+std::to_string(c));
        f.reference_asymmetry_ppb[c]=y(c,d+n)/y(c,n);
        if(!std::isfinite(f.reference_asymmetry_ppb[c])||f.reference_asymmetry_ppb[c]==0)throw std::runtime_error("Invalid reference asymmetry");
    }
    for(int i=0;i<n;i++){
        for(int c=0;c<components;c++){r[i]+=y(c,i);observed[i]+=y(c,d+i);f.rates(i,c)=y(c,i);}
        if(!(r[i]>0))throw std::runtime_error("Empty detector row "+std::to_string(i)+"; collect more events or enable missing detectors");
        observed[i]/=r[i];
        for(int c=0;c<components;c++){f.F(i,c)=y(c,i)/r[i];f.H(i,c)=y(c,d+i)/(r[i]*f.reference_asymmetry_ppb[c]);}
    }
    // Forecast covariance in ppb^2. It is distinct from finite-MC covariance V.
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        f.C(i,j)=b(i,j)/(polarization*polarization*days*86400*r[i]*r[j])*1e18;
        if(!std::isfinite(f.C(i,j)))throw std::runtime_error("Nonfinite projected covariance");
    }
    for(int i=0;i<n;i++){
        sd[i]=std::sqrt(f.C(i,i));if(!(sd[i]>0))throw std::runtime_error("Zero detector variance");
    }
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)f.correlation(i,j)=f.C(i,j)/(sd[i]*sd[j]);
    TDecompChol chol(f.correlation);
    if(!chol.Decompose())throw std::runtime_error("Singular detector covariance; insufficient statistics or redundant channels");
    const TMatrixD& u=chol.GetU(); // rho = U^T U; L = U^T.
    TMatrixD linv(n,n);
    for(int col=0;col<n;col++)for(int i=0;i<n;i++){
        double value=(i==col?1.:0.);for(int k=0;k<i;k++)value-=u(k,i)*linv(k,col);
        linv(i,col)=value/u(i,i);
    }
    TMatrixD weighted_h=f.H;
    for(int i=0;i<n;i++)for(int c=0;c<components;c++)weighted_h(i,c)/=sd[i];
    TMatrixD a=linv*weighted_h,scaled=a;
    TVectorD scales(components);
    for(int c=0;c<components;c++){
        for(int i=0;i<n;i++)scales[c]+=a(i,c)*a(i,c);
        scales[c]=std::sqrt(scales[c]);if(!(scales[c]>0))throw std::runtime_error("Missing fitted component");
        for(int i=0;i<n;i++)scaled(i,c)/=scales[c];
    }
    TDecompSVD svd(scaled);
    if(!svd.Decompose())throw std::runtime_error("SVD failed");
    const TVectorD& singular=svd.GetSig();
    if(singular[components-1]<=singular[0]*std::max(n,components)*std::numeric_limits<double>::epsilon())throw std::runtime_error("Rank deficient component fit");
    f.scaled_condition=singular[0]/singular[components-1];
    const TMatrixD& su=svd.GetU();const TMatrixD& sv=svd.GetV();
    TMatrixD pinv(components,n);
    for(int c=0;c<components;c++)for(int i=0;i<n;i++)for(int k=0;k<components;k++)pinv(c,i)+=sv(c,k)*su(i,k)/(singular[k]*scales[c]);
    TMatrixD g=pinv*linv;
    for(int c=0;c<components;c++)for(int i=0;i<n;i++)g(c,i)/=sd[i];
    for(int c=0;c<components;c++){
        double fitted=0;for(int i=0;i<n;i++)fitted+=g(c,i)*observed[i];
        if(std::abs(fitted-f.reference_asymmetry_ppb[c])>1e-8*(1+std::abs(f.reference_asymmetry_ppb[c])))throw std::runtime_error("Reference closure failed");
    }
    f.M.ResizeTo(components,components);f.M=transpose(a)*a;
    f.experimental_covariance.ResizeTo(components,components);f.experimental_covariance=pinv*transpose(pinv);
    // Propagate correlated dilution, asymmetry-shape and reference errors in H.
    TMatrixD j(n,components*dim);
    for(int c=0;c<components;c++)for(int i=0;i<n;i++){
        j(i,c*dim+i)=-observed[i]/r[i];j(i,c*dim+d+i)=1/r[i];
        j(i,c*dim+n)=f.H(i,c)*f.reference_asymmetry_ppb[c]/y(c,n);
        j(i,c*dim+d+n)=-f.H(i,c)/y(c,n);
    }
    TMatrixD derivative=g*j;
    f.mc_covariance.ResizeTo(components,components);f.mc_covariance=derivative*v*transpose(derivative);
    double max_mc=1;
    for(int i=0;i<components;i++)for(int k=i;k<components;k++){
        double mean=.5*(f.mc_covariance(i,k)+f.mc_covariance(k,i));f.mc_covariance(i,k)=f.mc_covariance(k,i)=mean;max_mc=std::max(max_mc,std::abs(mean));
    }
    for(int c=0;c<components;c++)if(f.mc_covariance(c,c)<-1e-9*max_mc)throw std::runtime_error("Negative propagated MC variance");
    f.combined_covariance.ResizeTo(components,components);f.combined_covariance=f.experimental_covariance+f.mc_covariance;
    f.F_mc_sigma.ResizeTo(n,components);
    for(int i=0;i<n;i++)for(int c=0;c<components;c++){
        double variance=0;
        for(int k=0;k<components;k++)for(int l=0;l<components;l++)variance+=((c==k?1.:0.)-f.F(i,c))*v(k*dim+i,l*dim+i)*((c==l?1.:0.)-f.F(i,c))/(r[i]*r[i]);
        f.F_mc_sigma(i,c)=std::sqrt(std::max(0.,variance));
    }
    for(const auto* matrix:{&f.F,&f.H,&f.C,&f.correlation,&f.M,&f.experimental_covariance,&f.mc_covariance,&f.combined_covariance,&f.F_mc_sigma})
        for(int i=0;i<matrix->GetNrows();i++)for(int j=0;j<matrix->GetNcols();j++)
            if(!std::isfinite((*matrix)(i,j)))throw std::runtime_error("Nonfinite fit output");
    return f;
}
}
