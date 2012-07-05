//nog enkele definities:
#ifdef PQ

#define __Q_CON

#endif

#ifdef PQG

#define __Q_CON
#define __G_CON

#endif

#ifdef PQGT1

#define __Q_CON
#define __G_CON
#define __T1_CON

#endif

#ifdef PQGT2

#define __Q_CON
#define __G_CON
#define __T2_CON

#endif

#ifdef PQGT

#define __Q_CON
#define __G_CON
#define __T1_CON
#define __T2_CON

#endif

#include "lapack.h"
#include "Matrix.h"
#include "BlockMatrix.h"
#include "Vector.h"
#include "BlockVector.h"

#include "LibInt.h"

#include "R.h"
#include "Gauss.h"
#include "input.h"

#include "preamble.h"

#include "CI_SPM.h"
#include "CI_TPM.h"
#include "CI_SPPM.h"
#include "CI_SPPM_m.h"
#include "CI_TPPM.h"

#include "CartInt.h"

#include "Transform.h"

#include "SI_SPM.h"
#include "SI_TPM.h"

#include "SphInt.h"

#include "Transform.h"

#include "Tools.h"

#include "SPM.h"
#include "TPM.h"
#include "PHM.h"
#include "DPM.h"
#include "PPHM.h"



#include "SUP.h"
#include "EIG.h"

#define HDF5_STATUS_CHECK(status) if(status < 0) std::cerr << __FILE__ << ":" << __LINE__ << ": Problem with writing to file. Status code=" << status << std::endl;

/* vim: set ts=3 sw=3 expandtab :*/
