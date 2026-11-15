#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <limits>

namespace {
bool armed=false,pending=false;
void require(bool condition,const char *message) {
    if(!condition) {std::fprintf(stderr,"GENESIS nonbond mirror: %s\n",message);std::abort();}
}
bool selected(const char *name,const char *expected) {
    const char *value=std::getenv(name);return value&&std::strcmp(value,expected)==0;
}
int mode() {
    static const int result=[] {
        const char *text=std::getenv("GENESIS_NONBOND_MIRROR");if(!text)return 0;
        char *end=nullptr;const long value=std::strtol(text,&end,10);
        require(end!=text&&!*end&&value>=0&&value<=2,"mode must be 0, 1 or 2");return int(value);
    }();return result;
}
}
extern "C" void genesis_nonbond_mirror_request(int eligible) {
    require(!armed&&!pending,"unconsumed prior force");
    // Only CHARMM short steps request retirement; full steps keep the CPU virial reader.
    armed=mode()&&eligible&&selected("GENESIS_FORCE_ASSEMBLY","2")&&selected("GENESIS_FORCE_BONDED","1")
        &&!selected("GENESIS_STEP_VERIFY","1")&&!selected("GENESIS_STEP_VERIFY","2");
}
extern "C" int genesis_nonbond_mirror_skip(void *host,std::size_t bytes,int virial,int real_bytes) {
    require(!pending,"prior GPU force was not assembled");
    const bool requested=armed;armed=false;
    if(!requested)return 0;
    require(!virial&&real_bytes==sizeof(double)&&host&&bytes&&bytes%sizeof(double)==0,"invalid short force request");
    if(selected("GENESIS_FORCE_VERIFY","1"))return 0;
    pending=true;
    if(mode()==2) {
        double *values=static_cast<double*>(host);
        for(std::size_t i=0;i<bytes/sizeof(double);++i)values[i]=std::numeric_limits<double>::quiet_NaN();
        std::printf("GENESIS_NONBOND_MIRROR: retired short force poisoned\n");
    }
    return 1;
}
extern "C" void genesis_nonbond_mirror_consume() {
    if(!pending)return;
    pending=false;
    if(mode()==2)std::printf("GENESIS_NONBOND_MIRROR: direct GPU force consumer pass\n");
}
