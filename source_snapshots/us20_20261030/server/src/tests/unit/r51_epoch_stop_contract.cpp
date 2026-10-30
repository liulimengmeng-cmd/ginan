#include "common/zhangEpochStopRequest.hpp"
#include <iostream>
#include <stdexcept>
#include <string>

static ZhangEpochStopRequest stop;
static void onSignal(int) {stop.request();}
static void check(bool result,const char* reason)
{if(!result)throw std::runtime_error(reason);}

int main()
{
    const int maxEpochs=-1;
    auto provenance=[&] {return std::string("effective_max_epochs=")+std::to_string(maxEpochs);};
    const auto before=provenance();
    check(!stop.pending() && stop.exitCode(false)==0,"clean run incorrectly interrupted");
    check(!stop.interruptedCaptureSchedule(true,true,false),"missing captures silently waived without a request");
    auto previous=std::signal(SIGINT,onSignal);
    std::raise(SIGINT);
    std::signal(SIGINT,previous);
    check(stop.pending() && before==provenance() && maxEpochs==-1,
          "signal mutated the authenticated configuration");
    check(stop.exitCode(false)==128+SIGINT,"interrupted run incorrectly marked complete");
    check(stop.interruptedCaptureSchedule(true,true,false),"intentional partial capture not identified");
    check(!stop.interruptedCaptureSchedule(true,true,true) && stop.exitCode(true)==EXIT_FAILURE,
          "signal waived a genuine checkpoint failure");
    check(!stop.interruptedCaptureSchedule(false,true,false) &&
          !stop.interruptedCaptureSchedule(true,false,false),"incorrect capture interruption status");
    check(!stop.request(),"second stop request not detected");
    std::cout<<"PASS actual SIGINT dispatch, immutable configuration, explicit interrupted exit, pending capture distinction, fatal precedence, repeated request\n";
}
