#include <stdio.h>
#include <nvector/nvector_serial.h>

#include "naunet.h"
#include "naunet_userdata.h"
#include "naunet_constants.h"
#include "fex.h"
#include "jtv.h"

int main() {

    realtype spy = 86400.0*365.0;

    UserData *data = new UserData();
    data->Tgas = 15.0;

    Naunet naunet;
    naunet.initSolver();

    realtype y[NSPECIES] = {0.0};
    y[0] = 0.4;
    y[1] = 0.4;
    y[2] = 0.1;
    y[3] = 0.1;

    naunet.solve(y, 10*spy, data);
    for (int i=0; i<NSPECIES; i++) {
        printf("\t%10.3e\n", y[i]);
    }

}

