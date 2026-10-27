// Real synthetic ROOT metadata tests. Scratch belongs in the caller's output area.
#include "../Deconvolution.C"
void test_root_inputs(const char* scratch){
    using namespace root_deconv;
    fs::path dir(scratch);fs::create_directories(dir);
    auto assign=[](remollTextFile& t,const std::string& name,const std::string& content){
        t.fFilenameSize=name.size()+1;t.fFilename=new char[t.fFilenameSize];std::strcpy(t.fFilename,name.c_str());
        t.fBufferSize=content.size()+1;t.fBuffer=new char[t.fBufferSize];std::strcpy(t.fBuffer,content.c_str());
    };
    auto make=[&](const char* filename,const std::string& mode,const std::string& phi,int seed,const char* fieldHash="",const char* extraSD=""){
        TFile file((dir/filename).string().c_str(),"RECREATE");TTree tree("T","Synthetic metadata test");
        double rate=1;tree.Branch("rate",&rate);tree.Fill();tree.Fill();tree.Write();
        remollRunData data;data.fNthrown=2;data.fSeed=seed;data.fGitInfo="synthetic-test-source";
        std::string macro="/remoll/evgen/set moller\n/remoll/beamcurr 65 uA\n/remoll/evgen/phmax "+phi+"\n/remoll/SD/enable_all\n";
        fs::path here=fs::absolute(fs::path(__FILE__)).parent_path().parent_path();
        fs::path channels=here/"data/channels.json";
        std::string txt=text(channels);std::regex pattern("\\\"detector\\\"\\s*:\\s*([0-9]+)");
        for(auto it=std::sregex_iterator(txt.begin(),txt.end(),pattern);it!=std::sregex_iterator();++it){
            std::string id=(*it)[1];
            if(!mode.empty())macro+="/remoll/SD/detect "+mode+" "+id+"\n";
            macro+="/remoll/SD/detect all "+id+"\n";
        }
        if(*fieldHash){filedata_t mag{};std::strcpy(mag.filename,"map.txt");std::strcpy(mag.hashsum,fieldHash);data.fMagData.push_back(mag);macro+="/remoll/addfield map.txt\n";}
        assign(data.fMacro,"test.mac",macro+extraSD+"/run/beamOn 2\n");
        remollTextFile geom;assign(geom,"geometry/test.gdml","<gdml/>");data.fGDMLFiles.push_back(geom);
        data.Write("run_data");file.Close();
    };
    make("one.root","boundaryhits","180 deg",1);make("two.root","surfacehits","360 deg",2);make("all.root","","180 deg",3);
    auto one=inspect(dir/"one.root",{110110},false,{});
    auto two=inspect(dir/"two.root",{110110},false,{});
    bool rejected=false;try{inspect(dir/"all.root",{110110},false,{});}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("detect all incorrectly establishes entry recording");
    make("disabled-range.root","boundaryhits","180 deg",6,"","/remoll/SD/disable_range 110000 112000\n");
    rejected=false;try{inspect(dir/"disabled-range.root",{110110},false,{});}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Disabled detector range accepted");
    make("enabled-range.root","","180 deg",7,"","/remoll/SD/disable_all\n/remoll/SD/enable_range 110110 110110\n/remoll/SD/detect_range BOUNDARYHITS 110110 110110\n");
    inspect(dir/"enabled-range.root",{110110},false,{});
    make("disabled-last.root","boundaryhits","180 deg",8,"","/remoll/SD/enable_range 110000 112000\n/remoll/SD/disable 110110\n");
    rejected=false;try{inspect(dir/"disabled-last.root",{110110},false,{});}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Detector command order ignored");
    for(auto key:{"phmin","phmax","thmin","thmax","thcommin","thcommax","emin","emax"}){
        Input changed=one;changed.settings[std::string("/remoll/evgen/")+key]="different";
        rejected=false;try{check_generation_settings(one,changed);}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("Mixed generation settings accepted");
    }
    check_generation_settings(one,one);Input different=two;different.process=1;check_generation_settings(one,different);
    // Cross-process samples must share transport, not just geometry and source.
    remollRunData fieldRun;fieldRun.fGitInfo="test-source";
    remollTextFile fieldGeometry;assign(fieldGeometry,"test.gdml","<gdml/>");fieldRun.fGDMLFiles.push_back(fieldGeometry);
    filedata_t field{};std::strcpy(field.filename,"/machine-a/map.txt");std::strcpy(field.hashsum,"same-content");
    fieldRun.fMagData.push_back(field);
    const std::string transportMacro="/remoll/addfield /machine-a/map.txt\n/remoll/field/scale map.txt 1\n/run/beamOn 2\n";
    auto original=transport_identity(fieldRun,transportMacro);
    std::strcpy(fieldRun.fMagData[0].filename,"/machine-b/map.txt");
    if(transport_identity(fieldRun,transportMacro)!=original)throw std::runtime_error("Field paths incorrectly change transport identity");
    std::strcpy(fieldRun.fMagData[0].hashsum,"changed-content");
    if(transport_identity(fieldRun,transportMacro)==original)throw std::runtime_error("Changed field contents accepted");
    std::strcpy(fieldRun.fMagData[0].hashsum,"same-content");
    if(transport_identity(fieldRun,"/remoll/rasx 2 mm\n"+transportMacro)==original)throw std::runtime_error("Changed beam settings accepted");
    if(transport_identity(fieldRun,"/remoll/field/scale map.txt 2\n"+transportMacro)==original)throw std::runtime_error("Ordered field commands lost");
    if(transport_identity(fieldRun,"/remoll/seed 99\n/remoll/beamcurr 1 uA\n/remoll/evgen/set elastic\n"+transportMacro)!=original)throw std::runtime_error("Generator/current/seed incorrectly change transport identity");
    different=one;different.transport="different";rejected=false;
    try{check_transport(one,different);}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Incompatible transport accepted");
    // ROOT-only users bypass Python's response validation: enforce it here too.
    const char* species[]={"e-","e+","gamma","mu-","mu+","pi-","pi+","neutron"};
    fs::path tables=dir/"tables";fs::create_directories(tables);
    auto writeTable=[&](const char* name,const std::string& rows){std::ofstream f(tables/(std::string("fit_param_xy_")+name+"_ifarm.csv"));f<<"energy,a,b,c,d,e\n"<<rows;};
    const std::string valid="1000,1,0,1,0,0\n1000,1,0,1,0,0\n2000,2,0,2,0,0\n";
    for(auto name:species)writeTable(name,valid);
    auto fit=portable_deconv::load(tables.string());
    if(fit.at(11).size()!=2)throw std::runtime_error("Identical upstream duplicate not handled");
    writeTable("pi-","1000,3,0,3,0,0\n");writeTable("pi+","1000,7,0,7,0,0\n");
    fit=portable_deconv::load(tables.string());
    if(fit.at(-211).at(1000)[0]!=3||fit.at(211).at(1000)[0]!=7)throw std::runtime_error("Charged pion PE tables interchanged");
    for(const std::string bad:{"1000,9,0,1,0,0\n","1500,broken,row\n","1500,nan,0,1,0,0\n","1500,1junk,0,1,0,0\n"}){
        writeTable("e-",valid+bad);rejected=false;
        try{portable_deconv::load(tables.string());}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("Malformed PE response accepted");
    }
    fs::path mainTables=dir/"main-tables",source=fs::absolute(fs::path(__FILE__)).parent_path().parent_path()/"data/main_pe";
    fs::create_directories(mainTables);for(const auto& f:fs::directory_iterator(source))fs::copy_file(f.path(),mainTables/f.path().filename(),fs::copy_options::overwrite_existing);
    MainModel validMain(mainTables);
    auto positionFile=mainTables/"r1_NewLookupTable.csv";const auto originalPosition=text(positionFile);
    auto quotedPosition=originalPosition;auto cellStart=quotedPosition.find('\n')+1;
    for(int i=0;i<8;i++)cellStart=quotedPosition.find(',',cellStart)+1;
    quotedPosition.insert(quotedPosition.find(',',cellStart),"\"");quotedPosition.insert(cellStart,"\"");
    {std::ofstream f(positionFile);f<<quotedPosition;}
    MainModel quotedMain(mainTables);
    if(quotedMain.position[0][0].p0!=validMain.position[0][0].p0)throw std::runtime_error("Quoted numeric PE coefficient changed value");
    {std::ofstream f(positionFile,std::ios::app);f<<"malformed,row\n";}rejected=false;
    try{MainModel bad(mainTables);}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Standalone main PE silently accepted malformed coefficients");
    {std::ofstream f(positionFile);f<<originalPosition;}
    auto energyFile=mainTables/"r1_EnergyDepTable.csv";auto energyText=text(energyFile);
    auto rowStart=energyText.find('\n')+1;energyText.replace(rowStart,energyText.find(',',rowStart)-rowStart,"9999");
    {std::ofstream f(energyFile);f<<energyText;}rejected=false;
    try{MainModel bad(mainTables);}catch(const std::exception&){rejected=true;}
    if(!rejected)throw std::runtime_error("Standalone main PE accepted unsupported energy grid");
    std::cout<<"PASS: transport identities and strict standalone ROOT PE tables\n";
    std::ofstream list(dir/"inputs.txt");list<<"one.root\ntwo.root\n";list.close();
    rejected=false;
    try{Deconvolution((dir/"inputs.txt").string().c_str(),(dir/"bad.root").string().c_str(),false);}
    catch(const std::exception& e){rejected=std::string(e.what()).find("Incompatible generator settings")!=std::string::npos;}
    if(!rejected||fs::exists(dir/"bad.root"))throw std::runtime_error("Public macro did not reject mixed domains before scanning");
    make("field0.root","boundaryhits","180 deg",4,"original-content");
    make("field1.root","boundaryhits","180 deg",5,"different-content");
    std::ofstream fieldList(dir/"fields.txt");fieldList<<"field0.root\nfield1.root\n";fieldList.close();rejected=false;
    try{Deconvolution((dir/"fields.txt").string().c_str(),(dir/"bad-field.root").string().c_str(),false);}
    catch(const std::exception& e){rejected=std::string(e.what()).find("Incompatible transport")!=std::string::npos;}
    if(!rejected||fs::exists(dir/"bad-field.root"))throw std::runtime_error("Public macro did not reject incompatible fields before scanning");
    std::cout<<"PASS: real ROOT metadata preserves entry recording; mixed generation domains rejected by public macro\n";
}
