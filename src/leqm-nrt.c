/*
   leqm-nrt is a  non-real-time implementation 
   of Leq(M) measurement according to ISO 21727:2004(E)
   "Cinematography -- Method of measurement of perceived
   loudness of motion-picture audio material" and the subsequent 
   revision ISO 21727:2016(E) "Cinematography — Method of measurement 
   of perceived loudness of short duration motion-picture audio material"

   Copyright (C) 2011-2013, 2017-2020 Luca Trisciani

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/




#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <iso646.h>

#ifdef _WIN32
#include <windows.h>
#elif defined __APPLE__
#include <sys/param.h>
#include <sys/sysctl.h>
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if  (defined(HAVE_LIBAVFORMAT)  &&  defined(HAVE_LIBAVCODEC) && defined(HAVE_LIBAVUTIL))
#define FFMPEG
#endif

#ifndef FFMPEG
#ifdef HAVE_LIBSNDFILE
#define SNDFILELIB
#else
#error no sndfile.h found
#endif
#endif

//#define DEBUG
//#define FFMPEG
//#define SNDFILELIB

#ifdef HAVE_LIBDI
#define DI
#else
#warning Compiling without Dolby Dialogue Intelligence
#endif


#ifdef FFMPEG
#include <stdio.h>
#include <libavutil/avassert.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/avutil.h>
//for calculation of true peak with ffmpeg I could use libavresample, see https://www.ffmpeg.org/doxygen/3.4/group__lavr.html
//see also ebur128.c for and example
#elif defined SNDFILELIB
#include <sndfile.h>
// Sample rate conversion src_process conflicts with Dolby DI
//#include <samplerate.h>
#include <stdint.h>
#endif


#ifdef DI
#include "di/di.h"

/* Input frame size, in samples */
//But this has nothing to do with the actual buffer used
// that is 2048ms with 75% overlap and di information
//actualization every 512ms
//The following is in samples, that is 400ms for 48kHz sampling
//#define INPUT_FRAME_SIZE            19200
// now I pass a variable buffersizesmpls to the function

/* Error codes */
#define ERR_INVALID_BIT_DEPTH       -1
#define ERR_CANNOT_OPEN_PCM_INPUT   -2
#define ERR_CANNOT_OPEN_OUTPUT      -3
#define ERR_CANNOT_OPEN_REFERENCE   -4
#define ERR_CANNOT_ALLOCATE_MEMORY  -5
#define ERR_INITIALIZATION_FAILED   -6

#endif


#define min(x, y) ({                \
		typeof(x) _min1 = (x);          \
		typeof(y) _min2 = (y);          \
		(void) (&_min1 == &_min2);      \
		_min1 < _min2 ? _min1 : _min2; })


#define max(x, y) ({                \
		typeof(x) _max1 = (x);          \
		typeof(y) _max2 = (y);          \
		(void) (&_max1 == &_max2);      \
		_max1 > _max2 ? _max1 : _max2; })


// When using DI there are some constraints as to the buffersize
// The same is true for LKFS

// Version 0.0.20 (C) Luca Trisciani 2011-2013, 2017-2020
// Tool from the DCP-Werkstatt Software Bundle



// COMPILATION
// compile for DEBUG with gcc -g -DEBUG -I/usr/include/ffmpeg -lfftw3 -lm -lpthread -lrt  -lavformat -lavcodec -lavutil -o leqm-nrt leqm-nrt.c

//otherwise  gcc -I/usr/include/ffmpeg -lm -lpthread -lrt -lavformat -lavcodec -lavutil -o leqm-nrt leqm-nrt.c


//this to do:
/*

   - with the github version and SNDFILELIB I get no memory leaks, so if I strip down of FFMPEG I should get the same here
   - It seams clean now also with FFMPEG see valgrind results

*/
/* Remember to avoid optimization if debugging, otherwise gdb behaviour could be unpredictable, jumping from line to line ...
   gcc -g3 -O0 -I/usr/include/di -lm -lpthread -lrt -lavformat -lavcodec -lavutil -L/usr/lib/di -o src/leqm-nrt  src/leqm-nrt.c -ldi -lrt -lpthread -lm -lavutil -lavformat -lavcodec
   */

#define MAX_RET 200000

//#define VERSION 20 //this is defined in config.h
/* //SRC conflicts with Dolby DI
#ifdef SNDFILELIB
SRC_DATA src_data;
#endif
*/

typedef struct
{

  double *vector;		//TruePeak in each channel
  int oversampling_ratio;
  int filtertaps;
  double *filter_coeffs;	//interpolation coefficients


} TruePeak;

struct Sum
{
  double csum;			// convolved sum
  double sum;			// flat sum
  int nsamples;
  double cmean;			//convolved mean
  double mean;
  double leqm;
  double rms;
  double lgleqm;		// level-gated leqm
  long int remainder_samples;	// stores number of samples in last non full buffer
#ifdef DI
  double dgleqm;		// dialogue-gated leqm
  double dialoguepercentual;	// Speech content %
#endif
};


typedef struct coefficient_context
{

  int sf;

  //HSF coefficient (K filter stage 1)

  double hsa0;
  double hsa1;
  double hsa2;
  double hsb0;
  double hsb1;
  double hsb2;

  //HPF coefficient (K filter stage 2 - RLB -> Revised Loudness B Weighting)

  double hpa0;
  double hpa1;
  double hpa2;
  double hpb0;
  double hpb1;
  double hpb2;

} coeff;



typedef struct LevelGate
{
  int gblocksize;		// in samples derived from sampling rate and block size in ms (no channel interleaving considered)
  int stepcounter;		//also keeps track of processed steps and total array length of results
  int totalsteps;		// at end of processing this will differ from stepcounter
  int shortperiods;		// totalsteps  should be equal to subdivs * shortperiods
  double ***LGresultarray;	//first index channels[i], second step short period[j] and third overlap step [k] (input for final calculation)
  int ops;			//offset per step in samples, derived from percentual overlap
  float overlap;		// percent
  int subdivs;			//this is necessary inside the worker to avoid cuncurrency
  float gatingThresholdFix;	// LKFS
  float relativeThresholdNegativeOffset;	// LKFS
  double *chgainconf;		// ITU 1770-4 gives only the 5.1 numbers
  int *chgateconf;		// for LKFS should be 3 for all channels except subwoofer
  double **oldBufferBackup;	//i is thread, j is sample - now this is still multichannel/interleaved
  //int  serialcheckstatus; // to serialize resource access
} LG;


typedef struct LevelGateLeqM
{
  int gblocksize;		// in samples derived from sampling rate and block size in ms (no channel interleaving considered)
  int stepcounter;		//also keeps track of processed steps and total array length of results
  int totalsteps;		// at end of processing this will differ from stepcounter
  int shortperiods;		// totalsteps  should be equal to subdivs * shortperiods
  double ***LGresultarray;	//first index channels[i], second step short period[j] and third overlap step [k] (input for final calculation)
  int ops;			//offset per step in samples, derived from percentual overlap
  float overlap;		// percent
  int subdivs;			//this is necessary inside the worker to avoid cuncurrency
  float gatingThresholdFix;	// LKFS
  float relativeThresholdNegativeOffset;	// LKFS
  double *chgainconf;		// ITU 1770-4 gives only the 5.1 numbers. Not really there is a long list
  int *chgateconf;		// for LKFS should be 3 for all channels except subwoofer
  double **oldBufferBackup;	//i is thread, j is sample
  //int serialcheckstatus; // to serialize resource access
  //int npcus; //number of cpus - 
} LGLeqM;




typedef struct Lg_Buffer
{
  int counter;
  double *bufferA;		//No Interleave 
  double *bufferB;
  double *bufferLG;
  double *bufferSwap;
} LG_Buf;



typedef struct Lg_BufferLeqM
{
  int counter;
  double *bufferA;		//No Interleave 
  double *bufferB;
  double *bufferLG;
  double *bufferSwap;
} LG_BufLeqM;


struct WorkerArgs
{
  double *argbuffer;
  float *di_argbuffer;
  int nsamples;
  int fullbuffer_nsamples;
  int nch;
  int npoints;
  double *ir;
  coeff *Kcoeffs;
  struct Sum *ptrtotsum;
  double *chconf;
  int *chgate;
  int shorttermindex;
  double **sc_shorttermarray;	// on the long run this should substitute shorttermarray. First index is channel and second is shorttermarray for a single channel 
  double *shorttermarray;	// here every shortperiod is the accumulation of all channels
#ifdef DI
  uint8_t **shorttermarray_di;	//this will point to the dialog / nodialog feature array shorttermdidecisionarray
  void **di_array;
#endif
  unsigned int *is_eof_array;
  unsigned int *buffered_samples_array;
  unsigned int *samples_read_array;
  int leqm10flag;
  int leqmlogflag;
  int lkfsflag;
  int leqmdiflag;
  int polyflag;
  int truepeakflag;
  TruePeak *truepeak;
  unsigned int sample_rate;	//needed by DI
  int channel;			//this is the channel being worked on at present. Needed by DI.
  LG *lg_ctx;
  LG_Buf *lg_buffers;
  LGLeqM *lg_ctx_leqmdi;
  LG_Buf *lg_buffers_leqmdi;
  int pthread_iteration;
  int ncpus;
  int worker_id;
#ifdef SNDFILELIB
  double src_output;
#endif
};



LG *LGCtx;
LGLeqM *LGCtxLeqMDI;

coeff *coeffs;

int precalculate_coeffs_K_filter (coeff * coeff_ctx, int samplerate);

int equalinterval (double *freqsamples, double *freqresp,
		   double *eqfreqsamples, double *eqfreqresp, int points,
		   int samplingfreq, int origpoints);
int equalinterval2 (double freqsamples[], double *freqresp,
		    double *eqfreqsamples, double *eqfreqresp, int points,
		    int samplingfreq, int origpoints, int bitdepthsoundfile);
int equalinterval3 (double freqsamples[], double *freqresp,
		    double *eqfreqsamples, double *eqfreqresp, int points,
		    int samplingfreq, int origpoints, int bitdepthsoundfile);
int convloglin (double *in, double *out, int points);
double convlinlog_single (double in);
double convloglin_single (double in);
int convolv_buff (double *sigin, double *sigout, double *impresp,
		  int sigin_dim, int impresp_dim);
double inputcalib (double dbdiffch);
int rectify (double *squared, double *inputsamples, int nsamples);
int accumulatech (double *chaccumulator, double *inputchannel, int nsamples);
double msaccumulate (double *inputbuffer, int nsamples);
#ifdef DI
int accumulatechwithdigate (double *chaccumulator, double *inputchannel,
			    int nsamples, int chgateconf,
			    struct WorkerArgs *workerargsinstance,
			    int channel_index);
#endif
int sumsamples (struct Sum *ts, double *inputsamples, double *cinputsamples,
		int nsamples);
int meanoverduration (struct Sum *oldsum);
void inversefft1 (double *eqfreqresp, double *ir, int npoints);
void inversefft2 (double *eqfreqresp, double *ir, int npoints);
void *worker_function (void *argstruct);
void *worker_function_gated2 (void *argstruct);
#ifdef DI
void *di_worker_function (void *argstruct);
#endif
void logleqm (FILE * filehandle, double featuretimesec, double temp_leqm);
double sumandshorttermavrg (double *channelaccumulator, int nsamples);
double logleqm10 (FILE * filehandle, double featuretimesec,
		  double longaverage);
#ifdef DI
void savedidecision (uint8_t ** dibytearray, int shortperiodidx,
		     uint8_t dibyte, int chnumb);
uint8_t **allocatedidecisionarray (int nchannels, int nshortperiods);
#endif
double **allocatesc_shorttermarray (int nchannels, int nshortperiods);
double ***allocateLGresultarray (int nchannels, int shortperiods,
				 int olsteps);
#ifdef DI
int freedidecisionarray (int nchannels,
			 uint8_t ** pt_shorttermdidecisionarray);
#endif
int freesc_shorttermarray (int nchannels, double **sc_shorttermarray);
int freeLGresultarray (int nchannels, int shortperiods, double ***lg_results);
unsigned int deintsamples (float *deintbuffer, int buffersizeframes,
			   double *intbuffer, int channel,
			   int tot_channel_numb);
#ifdef DI
int dolbydialogint (int channel, int tot_channels, uint8_t ** dibytearray,
		    double *buffer, unsigned int sample_rate,
		    int buffersizesmpls);
void print_di (int nchannels, int nshorttermperiods,
	       uint8_t ** di_decisionarray);
void dolbydifinalcomputation (double **sc_staa, uint8_t ** stdda, int nch,
			      int stpn, int chgateconfarray[],
			      struct Sum *ptSum);
void dolbydifinalcomputation2 (LGLeqM * pt_lgctx_leqmdi, int *pt_chgateconf,
			       int nchannels,
			       uint8_t ** stdda, double adlgthreshold);
#endif
double leqmtosum (double leqm, double ref);
void levelgatefinalcomputation (double **sc_staa, double linearthreshold, int nch, int stpn, struct Sum *ptSum);	//<-- Look at this 
void lkfs_finalcomputation (LG * pt_lgctx, int *pt_chgateconf, int nchannels);
void lkfs_finalcomputation_withdolbydi (LG * pt_lgctx, int *pt_chgateconf,
					int nchannels,
					uint8_t ** stdda,
					double adlgthreshold);

double *calc_lp_os_coeffs (int samplerate, int os_factor, int taps);
double truepeakcheck (double *in_buf, int ns, double truepeak, int os_ratio,
		      int filtertaps, double *coeff_vector);
TruePeak *init_truepeak_ctx (int ch, int os, int taps);
int freetruepeak (TruePeak * tp);

int calcSampleStepLG (float percentOverlap, int samplerate, int LGbufferms);
int K_filter_stage1 (double *smp_out, double *smp_in, int nsamples,
		     coeff * coeffctx);
int K_filter_stage2 (double *smp_out, double *smp_in, int nsamples,
		     coeff * coeffctx);
int M_filter (double *smp_out, double *smp_in, int samples, int samplerate);
LG_Buf *allocateLGBuffer (int samplenumber);
int freeLGBuffer (LG_Buf * pt_LG_Buf);

#ifdef FFMPEG

int transfer_decoded_data (AVFrame * ptr_frame,
			   struct WorkerArgs **dptrWorkerArgs, int w_id,
			   AVCodecContext * codecCon);
int transfer_decoded_samples (AVFrame * ptr_frame, double *buf,
			      AVCodecContext * codecCon, int buffersizs,
			      int *doublesample_index);
int transfer_remaining_decoded_samples (AVFrame * ptr_frame,
					double *bufremain,
					AVCodecContext * codecCon,
					int nxtsmpl, int *doublesample_index);
#endif

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t serialsignal = PTHREAD_COND_INITIALIZER;
pthread_cond_t serialsignal_leqmdi = PTHREAD_COND_INITIALIZER;

#ifdef DEBUG
void debuginit (double *pt_toarray, double value, int arraylength);
#endif



int
precalculate_coeffs_K_filter (coeff * coeff_ctx, int samplerate)
{

  /* For this calculation see (but the ITU standard has no indication and is older in the first 
     publication, so possibly there is also another source):

     "Algorithm described in Mansbridge, Stuart, Saoirse Finn, and Joshua D. Reiss. "Implementation
     and Evaluation of Autonomous Multi-track Fader Control." Paper presented at the
     132nd Audio Engineering Society Convention, Budapest, Hungary, 2012."
     (AES Convention Paper 8588) 

   */


  /* HSF prefilter coefficient calculation */

  double dBboost = 3.999843853973347;
  double cfreq = 1681.974450955533;
  double qf = 0.7071752369554196;

  double K = tan (M_PI * cfreq / samplerate);
  double Vh = pow (10.0, dBboost / 20.0);
  double Vb = pow (Vh, 0.4996667741545416);

  coeff_ctx->hsa0 = 1.0;
  double a_temp = 1.0 + K / qf + K * K;
  coeff_ctx->hsb0 = (Vh + Vb * K / qf + K * K) / a_temp;
  coeff_ctx->hsb1 = 2.0 * (K * K - Vh) / a_temp;
  coeff_ctx->hsb2 = (Vh - Vb * K / qf + K * K) / a_temp;
  coeff_ctx->hsa1 = 2.0 * (K * K - 1.0) / a_temp;
  coeff_ctx->hsa2 = (1.0 - K / qf + K * K) / a_temp;

  /* HPF calculations - alias Revised Loudness B Weighting */
  // reusing variables from preceding calculations

  cfreq = 38.13547087602444;
  qf = 0.5003270373238773;
  K = tan (M_PI * cfreq / samplerate);

  coeff_ctx->hpb0 = 1.0;
  coeff_ctx->hpb1 = -2.0;
  coeff_ctx->hpb2 = 1.0;
  coeff_ctx->hpa0 = 1.0;
  coeff_ctx->hpa1 = 2.0 * (K * K - 1.0) / (1.0 + K / qf + K * K);
  coeff_ctx->hpa2 = (1.0 - K / qf + K * K) / (1.0 + K / qf + K * K);

  return 0;


}












#ifdef FFMPEG








static double
get (uint8_t * a[], int ch, int index, int ch_count, enum AVSampleFormat f)
{
  const uint8_t *p;
  if (av_sample_fmt_is_planar (f))
    {
      f = av_get_alt_sample_fmt (f, 0);
      p = a[ch];
    }
  else
    {
      p = a[0];
      index = ch + index * ch_count;
    }

  switch (f)
    {
    case AV_SAMPLE_FMT_U8:
      return ((const uint8_t *) p)[index] / 127.0 - 1.0;
    case AV_SAMPLE_FMT_S16:
      return ((const int16_t *) p)[index] / 32767.0;
    case AV_SAMPLE_FMT_S32:
      return ((const int32_t *) p)[index] / 2147483647.0;
    case AV_SAMPLE_FMT_FLT:
      return ((const float *) p)[index];
    case AV_SAMPLE_FMT_DBL:
      return ((const double *) p)[index];
    default:
      av_assert0 (0);
    }
}


void
printAudioFrameInfo (const AVCodecContext * codecContext,
		     const AVFrame * frame)
{
  // See the following to know what data type (unsigned char, short, float, etc) to use to access the audio data:
  // http://ffmpeg.org/doxygen/trunk/samplefmt_8h.html#af9a51ca15301871723577c730b5865c5
  printf ("Audio frame info:\n");
  printf ("  Sample count: %d\n", frame->nb_samples);
  printf ("  Channel count: %d\n", codecContext->channels);
  printf ("  Format: %s\n",
	  av_get_sample_fmt_name (codecContext->sample_fmt));
  printf ("  Bytes per sample: %d\n",
	  av_get_bytes_per_sample (codecContext->sample_fmt));
  printf ("  Is planar? %d\n",
	  av_sample_fmt_is_planar (codecContext->sample_fmt));
  printf ("frame->linesize[0] tells you the size (in bytes) of each plane\n");

  if (codecContext->channels > AV_NUM_DATA_POINTERS
      && av_sample_fmt_is_planar (codecContext->sample_fmt))
    {
      printf
	("The audio stream (and its frames) have too many channels to fit in\n"
	 "frame->data. Therefore, to access the audio data, you need to use\n"
	 "frame->extended_data to access the audio data. It's planar, so\n"
	 "each channel is in a different element. That is:\n"
	 "  frame->extended_data[0] has the data for channel 1\n"
	 "  frame->extended_data[1] has the data for channel 2\n" "  etc.\n");
    }
  else
    {
      printf
	("Either the audio data is not planar, or there is enough room in\n"
	 "frame->data to store all the channels, so you can either use\n"
	 "frame->data or frame->extended_data to access the audio data (they\n"
	 "should just point to the same data).\n");
    }

  printf ("If the frame is planar, each channel is in a different element.\n"
	  "That is:\n"
	  "  frame->data[0]/frame->extended_data[0] has the data for channel 1\n"
	  "  frame->data[1]/frame->extended_data[1] has the data for channel 2\n"
	  "  etc.\n");

  printf ("If the frame is packed (not planar), then all the data is in\n"
	  "frame->data[0]/frame->extended_data[0] (kind of like how some\n"
	  "image formats have RGB pixels packed together, rather than storing\n"
	  " the red, green, and blue channels separately in different arrays.\n");
}
#endif

int
checkargvalue (const char *stringarg)
{
  if (stringarg == NULL)
    {
      printf ("Please provide required value after argument switch!\n");
      return 1;
    }
  else if (!(isdigit (stringarg[0])))
    {
      if ((strncmp (stringarg, "-", 1) == 0) && (isdigit (stringarg[1])))
	{
	  return 0;
	}
      else
	{
	  return 1;
	}
    }
  else
    {
      return 0;
    }
}


int
main (int argc, const char **argv)
{
  int npoints = 64;		// This value is low for precision. Calibration is done with 32768 points.
  int convpointsset = 0;
  int origpoints = 21;		//number of points in the standard CCIR filter
  int samplingfreq;		// this and the next is defined later taking it from sound file
  int bitdepth;
  // double normalizer;
  int timing = 0;
  int dolbydi = 0;
  int dolbydialt = 0;
  int levelgated = 0;
  int printdiinfo = 0;
  //double agsthreshold = 33.0; //Think about a sensible percentage value 
  /* the following variables are related to dolbydi */
  //int is_eof = 0;
  //int buffered_samples = 0;
  //int samples_read = 0;
  /* til here */
  struct timespec starttime;
  int fileopenstate = 0;
  int leqm10 = 0;
  int leqmlog = 0;
#if defined __unix__ || defined  __APPLE__
  int numCPU = sysconf (_SC_NPROCESSORS_ONLN) - 1;
#elif defined _WIN64 || defined _WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo (&sysinfo);
  int numCPU = sysinfo.dwNumberOfProcessors - 1;
#endif


  //ISO 21727:2004(E)
  // M Weighting
  double freqsamples[] =
    { 31, 63, 100, 200, 400, 800, 1000, 2000, 3150, 4000, 5000, 6300, 7100,
    8000, 9000, 10000, 12500, 14000, 16000, 20000, 31500
  };
  double freqresp_db[] =
    { -35.5, -29.5, -25.4, -19.4, -13.4, -7.5, -5.6, 0.0, 3.4, 4.9, 6.1, 6.6,
    6.4, 5.8, 4.5, 2.5, -5.6, -10.9, -17.3, -27.8, -48.3
  };


  double *eqfreqresp_db;
  eqfreqresp_db = NULL;
  double *eqfreqsamples;
  eqfreqsamples = NULL;
  double *eqfreqresp;
  eqfreqresp = NULL;
  double *ir;
  ir = NULL;

  double *channelconfcalvector;
  channelconfcalvector = NULL;
  int *channelgateconfvector;
  channelgateconfvector = NULL;
  printf
    ("leqm-nrt  Copyright (C) 2011-2013, 2017-2020 Luca Trisciani\nThis program comes with ABSOLUTELY NO WARRANTY,\nfor details on command line parameters see --help\nFirst argument is the audio file to be measured.\nOther parameters can follow in free order.\nThis is free software, and you are welcome to redistribute it\nunder the GPL v3 licence.\nProgram will use 1 + %d slave threads.\n",
     numCPU);
  //SndfileHandle file;
  //int oversampling_factor = 4;
#ifdef SNDFILELIB
  SNDFILE *file;
  file = NULL;
  SF_INFO sfinfo;
  memset (&sfinfo, 0, sizeof (sfinfo));
  int src_ratio = 8;		//fix at present, could become a parameter
#elif defined FFMPEG
  //av_register_all (); //it seems this is no more need for FFMPEG > 4.0
  AVFrame *frame = NULL;
  frame = av_frame_alloc ();
  AVFormatContext *formatContext = NULL;
  AVCodec *cdc = NULL;
  AVStream *audioStream = NULL;
  AVCodecContext *codecContext = NULL;


#endif

  FILE *leqm10logfile;
  leqm10logfile = NULL;
  FILE *leqmlogfile;
  leqmlogfile = NULL;
  int buffersizems = 850;	//ISO 21727:2004 do not contain any indication, TASA seems to indicate 1000, p. 8
  int buffersizesamples;	//samples per channel * channel * seconds
  int buffersizesamplesdi;	//samples per channel * seconds. This corresponds to INPUT_FRAME_SIZE
  double tempchcal[128];
  int tempchgate[128];
  tempchgate[0] = -1;		//this is just to see if the array is set
  int numcalread = 0;
  int numgateconfread = 0;
  double longperiod = 10.0;	//this is in minutes and correspond to the period for leqm10
  double threshold = 80.0;	//this is the threshold for the Allen metric
  double agsthreshold = 33.00;	// this is the speech percentage threshold for adaptive gate selection 
  double levelgatedthreshold = 38.00;	// this is the threshold for level gating
  double **sc_shorttermaveragedarray;	//first index is channel, second shorttermaveragedarray for each channel
  double *shorttermaveragedarray;
  /*

     One byte for 8 channels. If more than 8  channels are present more than 1 byte must be allocated.
     From little endian:
     Bit 7   6   5   4   3   2   1   0
     L   R   C  Sw  Ls  Rs LRs RRs

     No I changed my mind. As there is not so much need for memory, I will simply use an array of arrays.
     The first index is the channel and the second the short period index

   */
  uint8_t **shorttermdidecisionarray;

  /* At present I do not think anymore that it is necessary to have preprocessing. I think I will do without, but I need an array with so many DI instances as channels, array that I can point to in the argument to the workers */

  sc_shorttermaveragedarray = NULL;
  shorttermaveragedarray = NULL;
  shorttermdidecisionarray = NULL;

  int numbershortperiods = 0;
#ifdef FFMPEG
  int realnumbershortperiods = 0;
#endif
  int parameterstate = 0;
  int leqnw = 0;
  int lkfs = 0;
  int poly = 1;			//default to polynomial filtering starting from v. 0.20
  int truepeak = 0;
  int oversamp_ratio = 4;
  TruePeak *truepeak_ctx;


  char soundfilename[2048];
  // This is a requirement of sndfile library, do not forget it.

  const char helptext[] =
    "Order of parameters after audio file is free.\nPossible parameters are:\n--convpoints <integer number> \tUse convolution with n points interpolation instead of polynomial filter.\n\t\t\t\tDefault is polynomial filter.\n--numcpus <integer number> \tNumber of slave threads to speed up operation.\n--timing \t\t\tFor benchmarking speed.\n--chconfcal <dB correction> <dB correction> <etc. so many times as channels>\n--leqnw\t\t\t\toutput leq with no weighting\n--logleqm10\t\t\t(will also print Allen metric as output)\n--lkfs\t\t\t\tSwitch LKFS ITU 1770-4 on.\n--dolbydi\t\t\tSwitch Dolby Dialogue Intelligence on\n--chgateconf <0|1|2>, 0 = no gate, 1 = level gate (in dB) and 2 = dialogue gate\n--agsthreshold <speech %%>\tFor Leq(M,DI) and LKFS(DI) default 33%%.\n--levelgate <Leq(M)>\t\tThis will force level gating and deactivate speech gating\n--threshold <Leq(M)>\t\tThreshold used for Allen metric (default 80)\n--longperiod <minutes>\t\tLong period for leqm10 (default 10)\n--logleqm\t\t\tLog Leq(M) from start every buffersize ms.\n--buffersize <milliseconds>\t\t\tSize of Buffer in milliseconds.\n--truepeak\t\t\tShow true peak value\n--oversampling <n>\t\tDefault: 4 times\n--printdiinfo\t\t\tShow detailed speech intelligence information\n\nUsing:\ngnuplot -e \"plot \\\"logfile.txt\\\" u 1:2; pause -1\"\nit is possible to directly plot the logged data\n";


  if (argc == 1)
    {
      //printf(helptext);
      printf ("Please indicate a sound file to be processed.\n");
      return 0;
    }

  for (int in = 1; in < argc;)
    {
      if ((!(strncmp (argv[in], "-", 1) == 0)) && (argv[in] != NULL))
	{

#ifdef SNDFILELIB
	  if (fileopenstate == 0)
	    {
	      if (!(file = sf_open (argv[in], SFM_READ, &sfinfo)))
		{
		  printf
		    ("Error while opening audio file, could not open  %s\n.",
		     argv[in]);
		  puts (sf_strerror (NULL));
		  return 1;
		}


	      strcpy (soundfilename, argv[in]);
	      fileopenstate = 1;
	      printf ("Opened file: %s\n", argv[in]);
	      printf ("Sample rate: %d\n", sfinfo.samplerate);
	      printf ("Channels: %d\n", sfinfo.channels);
	      printf ("Format: %d\n", sfinfo.format);
	      printf ("Frames: %d\n", (int) sfinfo.frames);
	      channelconfcalvector =
		malloc (sizeof (double) * sfinfo.channels);
#ifdef DI
	      if (dolbydi || lkfs)
		{
#else
	      if (lkfs)
		{
#endif

		  channelgateconfvector =
		    malloc (sizeof (int) * sfinfo.channels);
		}

	      in++;
	      continue;
	    }
	  else
	    {
	      free (channelconfcalvector);
#ifdef DI
	      if (dolbydi || lkfs)
		{
#else
	      if (lkfs)
		{
#endif

		  free (channelgateconfvector);
		  channelgateconfvector = NULL;
		}
	      channelconfcalvector = NULL;

	      return 0;
	    }


#elif defined FFMPEG
	  if (fileopenstate == 0)
	    {
	      if (avformat_open_input (&formatContext, argv[in], NULL, NULL)
		  != 0)
		{
		  //av_free(frame);
		  av_frame_free (&frame);
		  printf ("Error opening the file\n");
		  return 1;
		}
	      //fileopenstate = 1;

	      if (avformat_find_stream_info (formatContext, NULL) < 0)
		{
		  //av_free(frame);
		  av_frame_free (&frame);
		  avformat_close_input (&formatContext);
		  printf ("Error finding the stream info\n");
		  return 1;
		}

	      // Find the audio stream
	      //AVCodec* cdc = NULL;
	      int streamIndex =
		av_find_best_stream (formatContext, AVMEDIA_TYPE_AUDIO, -1,
				     -1, &cdc, 0);
	      if (streamIndex < 0)
		{
		  //av_free(frame);
		  av_frame_free (&frame);
		  avformat_close_input (&formatContext);
		  printf ("Could not find any audio stream in the file\n");
		  return 1;
		}

	      audioStream = formatContext->streams[streamIndex];
	      //
	      //AVCodec * pCodec;
	      cdc = avcodec_find_decoder (audioStream->codecpar->codec_id);
	      codecContext = avcodec_alloc_context3 (cdc);
	      avcodec_parameters_to_context (codecContext,
					     audioStream->codecpar);
	      //avcodec_open2(codecContext, pCodec,NULL);


	      //
	      //codecContext = audioStream->codec;
	      //codecContext->codec = cdc;


	      /*
	       *
	       *
	       //pCodec=avcodec_find_decoder(pFormatCtx->streams[audioStream]->codecpar->codec_id);
	       pCodecCtx=avcodec_alloc_context3(pCodec);
	       avcodec_parameters_to_context(pCodecCtx,
	       pFormatCtx->streams[audioStream]->codecpar);
	       avcodec_open2(pCodecCtx, pCodec,NULL);
	       *
	       *
	       *
	       *
	       */

	      if (avcodec_open2 (codecContext, cdc, NULL) != 0)
		{
		  //      av_free(frame);
		  av_frame_free (&frame);
		  avformat_close_input (&formatContext);
		  printf ("Couldn't open the context with the decoder\n");
		  return 1;
		}

	      strcpy (soundfilename, argv[in]);
	      fileopenstate = 1;
	      printf ("This stream has %d ", codecContext->channels);
	      printf (" channels and a sample rate of %d ",
		      codecContext->sample_rate);
	      printf (" Hz\n");
	      printf ("The data is in the format %s\n",
		      av_get_sample_fmt_name (codecContext->sample_fmt));

	      channelconfcalvector =
		malloc (sizeof (double) * codecContext->channels);
#ifdef DI
	      if (dolbydi || lkfs)
		{
#else
	      if (lkfs)
		{
#endif
		  channelgateconfvector =
		    malloc (sizeof (double) * codecContext->channels);
		}

	      in++;
	      continue;
	    }
	  else
	    {			/*if (fileopenstate == 0) */
	      free (channelconfcalvector);
	      channelconfcalvector = NULL;
#ifdef DI
	      if (dolbydi || lkfs)
		{
#else
	      if (lkfs)
		{
#endif
		  free (channelgateconfvector);
		  channelgateconfvector = NULL;
		}
	      printf ("Please specify an input audio file.\n");
	      return 0;
	    }


#endif
	  continue;
	}			/*(!(strncmp(argv[in], "-", 1) == 0 */

      else if ((strcmp (argv[in], "--help") == 0) && (fileopenstate == 0))
	{
	  //in++;
	  printf (helptext);
	  return 0;

	}
      else if ((strcmp (argv[in], "--version") == 0) && (fileopenstate == 0))
	{

	  printf ("This is leqm-nrt version %s.\n", VERSION);

	  return 0;

	}
      else if (fileopenstate == 0)
	{
	  printf
	    ("Please indicate audio file to be processed as first argument.\n");
	  return 0;
	}

      if (strcmp (argv[in], "--chconfcal") == 0)
	{
	  /* as the order of parameter is free I have to postpone 
	     the check for consistency with the number of channels.
	     So first create a temporary array, whose number of element will be checked after 
	     the parsing of the command line parameters is finished. 
	     The calibration will be expressed in dB on the command line and converted to multiplier 
	     here so that it can be stored as a factor in the channelconfcalvector.
	   */

	  in++;
	  for (;;)
	    {
	      if (in < argc)
		{
		  //if (!(strncmp(argv[in], "-", 1) == 0)) { //changed this to allow negative numbers
		  if (!(strncmp (argv[in], "-", 1) == 0)
		      || isdigit (argv[in][1]))
		    {
		      tempchcal[numcalread++] = atof (argv[in++]);
		    }
		  else
		    break;
		}
	      else
		break;

	    }			//for

	  continue;
	}

      if (strcmp (argv[in], "--convpoints") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  npoints = atoi (argv[in + 1]);
	  in += 2;
	  printf ("Convolution points sets to %d.\n", npoints);
	  convpointsset = 1;
	  continue;

	}
      if (strcmp (argv[in], "--numcpus") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  numCPU = atoi (argv[in + 1]);
	  in += 2;
	  printf
	    ("Number of threads manually set to %d. Default is number of cores in the system minus one.\n",
	     numCPU);
	  continue;

	}
      if (strcmp (argv[in], "--truepeak") == 0)
	{
	  truepeak = 1;
	  in++;
	  printf ("True Peak value will be shown.\n");
	  continue;

	}
      if (strcmp (argv[in], "--oversampling") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  oversamp_ratio = atoi (argv[in + 1]);
	  in += 2;
	  if (truepeak == 0)
	    truepeak = 1;
	  printf
	    ("Oversampling ratio for truepeak set to %d. Also true peak reporting switched on.\n",
	     oversamp_ratio);
	  continue;

	}
#ifdef DI
      if (strcmp (argv[in], "--agsthreshold") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  agsthreshold = atof (argv[in + 1]);
	  in += 2;
	  printf ("Adaptive gate selection threshold set to %0.2f %%.\n",
		  agsthreshold);
	  continue;

	}
#endif

      if (strcmp (argv[in], "--chgateconf") == 0)
	{
	  /* as the order of parameter is free I have to postpone 
	     the check for consistency with the number of channels.
	     So first create a temporary array, whose number of element will be checked after 
	     the parsing of the command line parameters is finished. 
	   */

	  in++;
	  for (;;)
	    {
	      if (in < argc)
		{
		  //if (!(strncmp(argv[in], "-", 1) == 0)) { //changed this to allow negative numbers
		  if (isdigit (argv[in][0])
		      && ((atoi (argv[in]) == 0) || (atoi (argv[in]) == 1)
			  || (atoi (argv[in]) == 2)))
		    {		// 0 = no gate, 1 = level gate and 2 = dialogue gate. Level gate threasholds in dB are given on another array 
		      tempchgate[numgateconfread++] = atoi (argv[in++]);
		    }
		  else
		    break;
		}
	      else
		break;

	    }			//for

	  continue;
	}

      if (strcmp (argv[in], "--timing") == 0)
	{
	  timing = 1;
	  in++;
	  printf ("Execution time will be measured.\n");
	  continue;

	}

#ifdef DI
      if (strcmp (argv[in], "--dolbydi") == 0)
	{
	  dolbydi = 1;
	  in++;
	  printf ("Gating with Dolby Dialogue Intelligence.\n");
	  continue;

	}

      if (strcmp (argv[in], "--printdiinfo") == 0)
	{
	  printdiinfo = 1;
	  in++;
	  printf
	    ("Detailed speech intelligence information will be shown.\n");
	  continue;

	}

#endif
      if (strcmp (argv[in], "--lkfs") == 0)
	{
	  lkfs = 1;
	  in++;
	  printf ("Show ITU 1770-4 LKFS result.\n");
	  continue;

	}
      /*
         if (strcmp (argv[in], "--leqmdi") == 0)
         {
         lkfs = 1;
         in++;
         printf ("Show Leq(M,DI) measurement.\n");
         continue;

         } */

      /*
         if (strcmp (argv[in], "--polynomial") == 0)
         {
         poly = 1;
         in++;
         printf
         ("Using polynomial filter instead of convolution for M weighting.\n");
         continue;

         } */
      if (strcmp (argv[in], "--version") == 0)
	{
	  in++;
	  printf ("This is leqm-nrt version %s.\n", VERSION);
	  continue;

	}
      if (strcmp (argv[in], "--help") == 0)
	{
	  in++;
	  printf (helptext);
	  continue;

	}

      if (strcmp (argv[in], "--logleqm10") == 0)
	{
	  leqm10 = 1;
	  in++;
	  printf ("Leq(M,10m) data will be logged to the file leqm10.txt\n");
	  continue;

	}
      if (strcmp (argv[in], "--logleqm") == 0)
	{
	  leqmlog = 1;
	  in++;
	  printf ("Leq(M) data will be logged to the file leqmlog.txt\n");
	  continue;

	}

      if (strcmp (argv[in], "--leqnw") == 0)
	{
	  leqnw = 1;
	  in++;
	  printf ("Leq(nW) - unweighted -  will be outputted.\n");
	  continue;

	}

      if (strcmp (argv[in], "--buffersize") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  buffersizems = atoi (argv[in + 1]);
	  in += 2;
	  printf ("Buffersize will be set to %d milliseconds.\n",
		  buffersizems);
	  continue;

	}

      if (strcmp (argv[in], "--agsthreshold") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  agsthreshold = atof (argv[in + 1]);
	  in += 2;
	  printf
	    ("Adaptive gate selection based on speech percentage:  %.4f %%.\n",
	     agsthreshold);
	  continue;
	}

      if (strcmp (argv[in], "--levelgate") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  levelgatedthreshold = atof (argv[in + 1]);
	  in += 2;
	  levelgated = 1;
	  printf
	    ("Level gating threshold. This will force switching off speech gating and even adaptive gate selection. Specified threshold is %.4f Leq(M)\n",
	     levelgatedthreshold);
	  continue;
	}
      if (strcmp (argv[in], "--threshold") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  threshold = atof (argv[in + 1]);
	  in += 2;
	  printf ("Threshold for Allen metric set to %f Leq(M).\n",
		  threshold);
	  continue;

	}
      if (strcmp (argv[in], "--longperiod") == 0)
	{
	  if (checkargvalue (argv[in + 1]))
	    return 1;
	  longperiod = atof (argv[in + 1]);
	  in += 2;
	  printf ("Longperiod for Leq(M)X set to %f minutes.\n", longperiod);
	  continue;

	}

      if (parameterstate == 0)
	{
	  printf ("The command line switch you typed is not valid: %s\n",
		  argv[in]);
	  return -1;
	  break;
	}


    }				/* for (int in=1; in < argc;) */
  // Open audio file

  //postprocessing parameters

#ifdef SNDFILELIB

  if (numcalread == sfinfo.channels)
    {
      for (int cind = 0; cind < sfinfo.channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (tempchcal[cind]);

	}
    }
  else if ((numcalread == 0) && (sfinfo.channels == 2))
    {
      double conf20[2] = { 0, 0 };
      for (int cind = 0; cind < sfinfo.channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf20[cind]);
	}
      printf
	("Using input channel calibration for 2.0 configuration:\n0 0\n");
    }
  else if ((numcalread == 0) && (sfinfo.channels == 6))
    {
      double conf51[6] = { 0, 0, 0, 0, -3, -3 };
      for (int cind = 0; cind < sfinfo.channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf51[cind]);
	}
      printf
	("Using input channel calibration for 5.1 configuration:\n0 0 0 0 -3 -3\n");
    }
  else if ((numcalread == 0) && (sfinfo.channels == 8))
    {
      double conf71[8] = { 0, 0, 0, 0, -3, -3, -3, -3 };
      for (int cind = 0; cind < sfinfo.channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf71[cind]);
	}
      printf
	("Using input channel calibration for 7.1 configuration:\n0 0 0 0 -3 -3 -3 -3\n");

    }

#elif defined FFMPEG

  if (numcalread == codecContext->channels)
    {
      for (int cind = 0; cind < codecContext->channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (tempchcal[cind]);
	}
    }

  else if ((numcalread == 0) && (codecContext->channels == 2))
    {
      double conf20[2] = { 0, 0 };
      for (int cind = 0; cind < codecContext->channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf20[cind]);
	}
      printf
	("Using input channel calibration for 2.0 configuration:\n0 0\n");
    }
  else if ((numcalread == 0) && (codecContext->channels == 6))
    {
      double conf51[6] = { 0, 0, 0, 0, -3, -3 };
      for (int cind = 0; cind < codecContext->channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf51[cind]);
	}
      printf
	("Using input channel calibration for 5.1 configuration:\n0 0 0 0 -3 -3\n");
    }
  else if ((numcalread == 0) && (codecContext->channels == 8))
    {
      double conf71[8] = { 0, 0, 0, 0, -3, -3, -3, -3 };
      for (int cind = 0; cind < codecContext->channels; cind++)
	{
	  channelconfcalvector[cind] = convloglin_single (conf71[cind]);
	}
      printf
	("Using input channel calibration for 7.1 configuration:\n0 0 0 0 -3 -3 -3 -3\n");
    }
#endif
  else
    {

      printf
	("Either you specified a different number of calibration than number of channels in the file or you do not indicate any calibration and the program cannot infer one from the number of channels. Please specify a channel calibration on the command line.\n");

      free (channelconfcalvector);
      channelconfcalvector = NULL;
      return 0;
    }

  if (truepeak)
    {
      int filtertaps = 12 * oversamp_ratio;
#ifdef SNDFILELIB
      switch (sfinfo.samplerate)
#elif defined FFMPEG
      switch (codecContext->sample_rate)
#endif
	{
	case 44100:
	  oversamp_ratio = oversamp_ratio;
	  filtertaps = filtertaps;
	  break;
	case 48000:
	  oversamp_ratio = oversamp_ratio;
	  filtertaps = filtertaps;
	  break;
	case 96000:
	  oversamp_ratio = oversamp_ratio / 2;
	  filtertaps = filtertaps / 2;
	  break;
	default:
	  oversamp_ratio = 1;	// no oversampling for theoretical higher sampling rate
	  filtertaps = filtertaps / 4;	// would this work? filter should not execute ! 
	  break;
	}


      //initialize oversampling filter
#ifdef FFMPEG
      truepeak_ctx =
	init_truepeak_ctx (codecContext->channels, oversamp_ratio,
			   filtertaps);
      truepeak_ctx->filter_coeffs =
	calc_lp_os_coeffs (codecContext->sample_rate, oversamp_ratio,
			   filtertaps);

#elif defined SNDFILELIB
      truepeak_ctx =
	init_truepeak_ctx (sfinfo.channels, oversamp_ratio, filtertaps);
      truepeak_ctx->filter_coeffs =
	calc_lp_os_coeffs (sfinfo.samplerate, oversamp_ratio, filtertaps);
#endif
    }



  if (leqm10)
    {
      char tempstring[2048];
      strcpy (tempstring, soundfilename);
      strcat (tempstring, ".leqm10.txt");
      leqm10logfile = fopen (tempstring, "w");
      if (leqm10logfile == NULL)
	{
	  printf ("Could not open file to write log leqm10 data!\n");
	}
    }

#ifdef DI
  if (dolbydi)
    {
      if (tempchgate[0] == -1)
	{
	  printf
	    ("Dolby Dialog Intelligence is active for channel L C R. Please\nuse --chgateconf <0|1|2> ... for other configurations.\n");
#ifdef SNDFILELIB
	  if ((sfinfo.channels == 6) || (sfinfo.channels < 3))
	    {
#elif defined FFMPEG
	  if ((codecContext->channels == 6) || (codecContext->channels < 3))
	    {
#endif
	      tempchgate[0] = 2;
	      tempchgate[1] = 2;
	      tempchgate[2] = 2;
	      tempchgate[3] = 0;
	      tempchgate[4] = 0;
	      tempchgate[5] = 0;

#ifdef SNDFILELIB
	    }
	  else if (sfinfo.channels == 8)
	    {
#elif defined FFMPEG
	    }
	  else if (codecContext->channels == 8)
	    {
#endif
	      tempchgate[0] = 2;
	      tempchgate[1] = 2;
	      tempchgate[2] = 2;
	      tempchgate[3] = 0;
	      tempchgate[4] = 0;
	      tempchgate[5] = 0;
	      tempchgate[6] = 0;
	      tempchgate[7] = 0;

	    }
	  else
	    {
	      printf
		("Please specify which channel(s) should go through Dialog Intelligence.\n");
	      return 0;
	    }
	}
    }


#endif


  if (leqmlog)
    {
      char tempstring[2048];
      strcpy (tempstring, soundfilename);
      strcat (tempstring, ".leqmlog.txt");
      leqmlogfile = fopen (tempstring, "w");
      if (leqmlogfile == NULL)
	{
	  printf ("Could not open file to write log leqm data!\n");
	}
    }


  if (timing)
    {

      clock_gettime (CLOCK_MONOTONIC, &starttime);
    }

  if (convpointsset)
    {
      poly = 0;
      printf ("Using convolution instead of polynomial filtering.\n");
    }

  if (lkfs)
    {
      printf
	("Setting buffersize to 400ms as per ITU 1770-4 requirements for LKFS measurement.\n");
      buffersizems = 400;
      //initialize coeff struct and calculates coefficients for K filter stage 1 and 2
      coeffs = malloc (sizeof (coeff));
#ifdef SNDFILELIB
      precalculate_coeffs_K_filter (coeffs, sfinfo.samplerate);
#elif defined FFMPEG
      precalculate_coeffs_K_filter (coeffs, codecContext->sample_rate);
#endif


    }



  // reading to a double or float buffer with sndfile take care of normalization
  /*
     static double  buffer[BUFFER_LEN]; // it seems this must be static. I don't know why
   */
  double *buffer;
#ifdef FFMPEG
  double *remainbuffer;
#endif
  // buffer = new double [BUFFER_LEN];
  //buffersizesamples = (sfinfo.samplerate*sfinfo.channels*buffersizems)/1000;

#ifdef SNDFILELIB
  if ((sfinfo.samplerate * buffersizems) % 1000)
    {
#elif defined FFMPEG
  if ((codecContext->sample_rate * buffersizems) % 1000)
    {
#endif
      printf
	("Please fine tune the buffersize according to the sampling rate\n");
      //close file
      // free memory
      // write a function to do that 
      return 1;
    }


#ifdef SNDFILELIB
  buffersizesamples =
    (sfinfo.samplerate * sfinfo.channels * buffersizems) / 1000;
  buffer = malloc (sizeof (double) * buffersizesamples);
  buffersizesamplesdi = (sfinfo.samplerate * buffersizems) / 1000;
  samplingfreq = sfinfo.samplerate;

#ifdef DI
  if (leqm10 || leqmlog || dolbydi)
    {
#else
  if (leqm10 || leqmlog)
    {
#endif

#ifdef DI
      //if duration < 10 mm exit
      if (leqm10)
	{
#endif
	  double featdursec = sfinfo.frames / sfinfo.samplerate;
	  if ((featdursec / 60.0) < longperiod)
	    {
	      printf ("The audio file is too short to measure Leq(m10).\n");
	      //return 0;
	    }
#ifdef DI
	}
#endif

      //how many short periods in overall duration //ATTENTION this calculation may return 0  for buffer sizes of less than 10 ms.
      //Maybe cast to double!!
      //Also remember that this is necessary also for dolbydi
      int remainderl =
	sfinfo.frames % (sfinfo.samplerate * buffersizems / 1000);
      if (remainderl == 0)
	numbershortperiods =
	  sfinfo.frames / (sfinfo.samplerate * buffersizems / 1000);
      else
	numbershortperiods =
	  sfinfo.frames / (sfinfo.samplerate * buffersizems / 1000) + 1;

      //allocate array // well the first is wrong

      shorttermaveragedarray =
	malloc (sizeof (*shorttermaveragedarray) * numbershortperiods);
#ifdef DI
      if (dolbydi)
	{
	  shorttermdidecisionarray = allocatedidecisionarray (sfinfo.channels, numbershortperiods);	/* Allocate DI array for all channels */
	  sc_shorttermaveragedarray =
	    allocatesc_shorttermarray (sfinfo.channels, numbershortperiods);
	}
#endif
    }				/* if (leqm10) */




#elif defined FFMPEG
  buffersizesamples =
    (codecContext->sample_rate * codecContext->channels * buffersizems) /
    1000;
  buffer = malloc (sizeof (double) * buffersizesamples);
  remainbuffer = malloc (sizeof (double) * buffersizesamples);
  buffersizesamplesdi = (codecContext->sample_rate * buffersizems) / 1000;
  samplingfreq = codecContext->sample_rate;
  // postpone this because I cannot get duration or number of frames
  // but I cannot postpone this!!
  //it seems I cannot get total number of audio frames in ffmpeg so I will simply allocate enough
  //memory for a 5 hours feature, ok?
#ifdef DI
  if (leqm10 || leqmlog || dolbydi)
    {
#else
  if (leqm10 || leqmlog)
    {
#endif
      //numbershortperiods = (int) (180000.00 / (((double) codecContext->sample_rate) * (double) buffersizems/1000.00) + 1); //this is wrong, because 180000 cannot be be number of frames, also why should we devide by the sample rate if 18000 is just seconds?
      numbershortperiods =
	(int) (18000.00 / ((double) buffersizems / 1000.00) + 1);
      shorttermaveragedarray =
	malloc (sizeof (*shorttermaveragedarray) * numbershortperiods);
#ifdef DI
      if (dolbydi)
	{
	  shorttermdidecisionarray = allocatedidecisionarray (codecContext->channels, numbershortperiods);	/* Allocate DI array for all channels */
	  sc_shorttermaveragedarray =
	    allocatesc_shorttermarray (codecContext->channels,
				       numbershortperiods);
	}
#endif
    }				/* if (leqm10) */

#endif


  if (lkfs)
    {

      LGCtx = malloc (sizeof (LG));	// this must be freed
      LGCtx->overlap = 75.0;	//pecentage overlap. should become a command line switch later, default to this value as ITU 1770-4.
      LGCtx->stepcounter = 0;
      LGCtx->gatingThresholdFix = -70.0;
      LGCtx->relativeThresholdNegativeOffset = -10.0;
      //LGCtx->serialcheckstatus = 0;

#ifdef SNDFILELIB
      //calculate step offset from percent overlap
      LGCtx->ops =
	calcSampleStepLG (LGCtx->overlap, sfinfo.samplerate, buffersizems);
      LGCtx->gblocksize = (sfinfo.samplerate * buffersizems / 1000);	//At present this is the same as buffersizems, but in samples (no interleaving)




      //first calculate or guestimate total step number, see stepcounter in LG 
      int remainder = sfinfo.frames % LGCtx->ops;

      /* total step calculation */

      if (remainder == 0)
	{
	  LGCtx->totalsteps = sfinfo.frames / LGCtx->ops;
	}
      else
	LGCtx->totalsteps = sfinfo.frames / LGCtx->ops + 1;

      remainder = LGCtx->gblocksize % LGCtx->ops;

      if (!(remainder))
	{
	  LGCtx->subdivs = LGCtx->gblocksize / LGCtx->ops;
	  //LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs;
	}
      else
	{
	  LGCtx->subdivs = LGCtx->gblocksize / LGCtx->ops + 1;
	  //  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs + 1;
	}

      remainder = LGCtx->totalsteps % LGCtx->subdivs;
      if (!remainder)
	{
	  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs;
	}
      else
	{

	  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs + 1;
	}

      LGCtx->LGresultarray = allocateLGresultarray (sfinfo.channels, LGCtx->shortperiods, LGCtx->subdivs);	//numbershortperiods is not the real length of the data, see above
      LGCtx->chgainconf = malloc (sizeof (double) * sfinfo.channels);
      LGCtx->chgateconf = malloc (sizeof (double) * sfinfo.channels);
      //maybe
      //other speaker configurations?
      if (sfinfo.channels == 6)
	{
	  LGCtx->chgainconf[0] = 1;
	  LGCtx->chgainconf[1] = 1;
	  LGCtx->chgainconf[2] = 1;
	  LGCtx->chgainconf[3] = 0;
	  LGCtx->chgainconf[4] = 1.41;
	  LGCtx->chgainconf[5] = 1.41;

	  LGCtx->chgateconf[0] = 3;
	  LGCtx->chgateconf[1] = 3;
	  LGCtx->chgateconf[2] = 3;
	  LGCtx->chgateconf[3] = 0;
	  LGCtx->chgateconf[4] = 3;
	  LGCtx->chgateconf[5] = 3;



	}
      else
	{
	  for (int i = 0; i < sfinfo.channels; i++)
	    {
	      LGCtx->chgainconf[i] = 1;
	      if (i != 3)
		{
		  LGCtx->chgateconf[i] = 3;
		}
	      else
		{
		  LGCtx->chgateconf[i] = 0;
		}
	    }
	}
      //first calculate or guestimate total step number, see stepcounter in LG
      //allocate overlap result array, see olresultarrey in LG



#elif defined FFMPEG

      //calculate step offset from percent overlap
      LGCtx->ops =
	calcSampleStepLG (LGCtx->overlap, codecContext->sample_rate,
			  buffersizems);
      LGCtx->gblocksize = (codecContext->sample_rate * buffersizems / 1000);	//At present this is the same as buffersizems, but in samples (no interleaving)
      //first calculate or guestimate total step number, see stepcounter in LG 
      LGCtx->totalsteps =
	((int)
	 (18000.00 /
	  (((double) LGCtx->ops) / ((double) codecContext->sample_rate)))) +
	1;

      int remainder = LGCtx->gblocksize % LGCtx->ops;
      /*subdivs */
      if (!(remainder))
	{
	  LGCtx->subdivs = LGCtx->gblocksize / LGCtx->ops;
	}
      else
	{
	  LGCtx->subdivs = LGCtx->gblocksize / LGCtx->ops + 1;
	}
      /* shortperiods */

      remainder = LGCtx->totalsteps % LGCtx->subdivs;
      if (!remainder)
	{
	  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs;
	}
      else
	{
	  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs + 1;

	}
      LGCtx->LGresultarray = allocateLGresultarray (codecContext->channels, LGCtx->shortperiods, LGCtx->subdivs);	//numbershortperiods is not the real length of data
      //first calculate or guestimate total step number, see stepcounter in LG
      //allocate overlap result array, see olresultarrey in LG

      LGCtx->chgainconf = malloc (sizeof (double) * codecContext->channels);
      LGCtx->chgateconf = malloc (sizeof (double) * codecContext->channels);
      //maybe 
      if (codecContext->channels == 6)
	{
	  LGCtx->chgainconf[0] = 1;
	  LGCtx->chgainconf[1] = 1;
	  LGCtx->chgainconf[2] = 1;
	  LGCtx->chgainconf[3] = 0;
	  LGCtx->chgainconf[4] = 1.41;
	  LGCtx->chgainconf[5] = 1.41;

	  LGCtx->chgateconf[0] = 3;
	  LGCtx->chgateconf[1] = 3;
	  LGCtx->chgateconf[2] = 3;
	  LGCtx->chgateconf[3] = 0;
	  LGCtx->chgateconf[4] = 3;
	  LGCtx->chgateconf[5] = 3;

	  //remember to costrain buffersizems to 400ms in case lkfs is used
	}
      else
	{
	  for (int i = 0; i < codecContext->channels; i++)
	    {
	      LGCtx->chgainconf[i] = 1;
	      if (i != 3)
		{
		  LGCtx->chgateconf[i] = 3;
		}
	      else
		{
		  LGCtx->chgateconf[i] = 0;
		}
	    }
	}
#endif

#ifdef SNDFILELIB
      //LGCtx->oldBufferBackup = malloc (sizeof (double *) * sfinfo.channels);
      LGCtx->oldBufferBackup = malloc (sizeof (double *) * numCPU);
      //for (int i = 0; i < sfinfo.channels; i++)
      for (int i = 0; i < numCPU; i++)
	{
	  LGCtx->oldBufferBackup[i] =
	    //malloc (sizeof (double) * LGCtx->gblocksize);
	    malloc (sizeof (double) * LGCtx->gblocksize * sfinfo.channels);
	}
#elif defined FFMPEG
      //LGCtx->oldBufferBackup =
      //  malloc (sizeof (double *) * codecContext->channels);
      LGCtx->oldBufferBackup = malloc (sizeof (double *) * numCPU);

      //    for (int i = 0; i < codecContext->channels; i++)
      for (int i = 0; i < numCPU; i++)
	{
	  *(LGCtx->oldBufferBackup + i) =
	    malloc (sizeof (double) * LGCtx->gblocksize *
		    codecContext->channels);
	}
#endif
    }				/* if (lkfs) */

#ifdef DI
  /* LEQMDI Insert */


  if (dolbydi)
    {

      LGCtxLeqMDI = malloc (sizeof (LGLeqM));	// this must be freed
      LGCtxLeqMDI->overlap = 75.0;	//pecentage overlap. should become a command line switch later, default to this value as ITU 1770-4.
      LGCtxLeqMDI->stepcounter = 0;
      LGCtxLeqMDI->gatingThresholdFix = -70.0;
      LGCtxLeqMDI->relativeThresholdNegativeOffset = -10.0;


#ifdef SNDFILELIB
      //calculate step offset from percent overlap
      LGCtxLeqMDI->ops =
	calcSampleStepLG (LGCtxLeqMDI->overlap, sfinfo.samplerate,
			  buffersizems);
      LGCtxLeqMDI->gblocksize = (sfinfo.samplerate * buffersizems / 1000);	//At present this is the same as buffersizems, but in samples (no interleaving)




      //first calculate or guestimate total step number, see stepcounter in LG 
      int remainder = sfinfo.frames % LGCtxLeqMDI->ops;

      /* total step calculation */

      if (remainder == 0)
	{
	  LGCtxLeqMDI->totalsteps = sfinfo.frames / LGCtxLeqMDI->ops;
	}
      else
	LGCtxLeqMDI->totalsteps = sfinfo.frames / LGCtxLeqMDI->ops + 1;

      remainder = LGCtxLeqMDI->gblocksize % LGCtxLeqMDI->ops;

      if (!(remainder))
	{
	  LGCtxLeqMDI->subdivs = LGCtxLeqMDI->gblocksize / LGCtxLeqMDI->ops;
	  //LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs;
	}
      else
	{
	  LGCtxLeqMDI->subdivs =
	    LGCtxLeqMDI->gblocksize / LGCtxLeqMDI->ops + 1;
	  //  LGCtx->shortperiods = LGCtx->totalsteps / LGCtx->subdivs + 1;
	}

      remainder = LGCtxLeqMDI->totalsteps % LGCtxLeqMDI->subdivs;
      if (!remainder)
	{
	  LGCtxLeqMDI->shortperiods =
	    LGCtxLeqMDI->totalsteps / LGCtxLeqMDI->subdivs;
	}
      else
	{

	  LGCtxLeqMDI->shortperiods =
	    LGCtxLeqMDI->totalsteps / LGCtxLeqMDI->subdivs + 1;
	}

      LGCtxLeqMDI->LGresultarray = allocateLGresultarray (sfinfo.channels, LGCtxLeqMDI->shortperiods, LGCtxLeqMDI->subdivs);	//numbershortperiods is not the real length of the data, see above
      LGCtxLeqMDI->chgainconf = malloc (sizeof (double) * sfinfo.channels);
      LGCtxLeqMDI->chgateconf = malloc (sizeof (double) * sfinfo.channels);
      //maybe
      //other speaker configurations?
      if (sfinfo.channels == 6)
	{
	  LGCtxLeqMDI->chgainconf[0] = 1;
	  LGCtxLeqMDI->chgainconf[1] = 1;
	  LGCtxLeqMDI->chgainconf[2] = 1;
	  LGCtxLeqMDI->chgainconf[3] = 0;
	  LGCtxLeqMDI->chgainconf[4] = 1;
	  LGCtxLeqMDI->chgainconf[5] = 1;

	  LGCtxLeqMDI->chgateconf[0] = 3;
	  LGCtxLeqMDI->chgateconf[1] = 3;
	  LGCtxLeqMDI->chgateconf[2] = 3;
	  LGCtxLeqMDI->chgateconf[3] = 0;
	  LGCtxLeqMDI->chgateconf[4] = 3;
	  LGCtxLeqMDI->chgateconf[5] = 3;



	}
      else
	{
	  for (int i = 0; i < sfinfo.channels; i++)
	    {
	      LGCtxLeqMDI->chgainconf[i] = 1;
	      if (i != 3)
		{
		  LGCtxLeqMDI->chgateconf[i] = 3;
		}
	      else
		{
		  LGCtxLeqMDI->chgateconf[i] = 0;
		}
	    }
	}
      //first calculate or guestimate total step number, see stepcounter in LG
      //allocate overlap result array, see olresultarrey in LG



#elif defined FFMPEG

      //calculate step offset from percent overlap
      LGCtxLeqMDI->ops =
	calcSampleStepLG (LGCtxLeqMDI->overlap, codecContext->sample_rate,
			  buffersizems);
      LGCtxLeqMDI->gblocksize = (codecContext->sample_rate * buffersizems / 1000);	//At present this is the same as buffersizems, but in samples (no interleaving)
      //first calculate or guestimate total step number, see stepcounter in LG 
      LGCtxLeqMDI->totalsteps =
	((int)
	 (18000.00 /
	  (((double) LGCtxLeqMDI->ops) /
	   ((double) codecContext->sample_rate)))) + 1;

      int remainder = LGCtxLeqMDI->gblocksize % LGCtxLeqMDI->ops;
      /*subdivs */
      if (!(remainder))
	{
	  LGCtxLeqMDI->subdivs = LGCtxLeqMDI->gblocksize / LGCtxLeqMDI->ops;
	}
      else
	{
	  LGCtxLeqMDI->subdivs =
	    LGCtxLeqMDI->gblocksize / LGCtxLeqMDI->ops + 1;
	}
      /* shortperiods */

      remainder = LGCtxLeqMDI->totalsteps % LGCtxLeqMDI->subdivs;
      if (!remainder)
	{
	  LGCtxLeqMDI->shortperiods =
	    LGCtxLeqMDI->totalsteps / LGCtxLeqMDI->subdivs;
	}
      else
	{
	  LGCtxLeqMDI->shortperiods =
	    LGCtxLeqMDI->totalsteps / LGCtxLeqMDI->subdivs + 1;

	}
      LGCtxLeqMDI->LGresultarray = allocateLGresultarray (codecContext->channels, LGCtxLeqMDI->shortperiods, LGCtxLeqMDI->subdivs);	//numbershortperiods is not the real length of data
      //first calculate or guestimate total step number, see stepcounter in LG
      //allocate overlap result array, see olresultarrey in LG

      LGCtxLeqMDI->chgainconf =
	malloc (sizeof (double) * codecContext->channels);
      LGCtxLeqMDI->chgateconf =
	malloc (sizeof (double) * codecContext->channels);
      //maybe 
      if (codecContext->channels == 6)
	{
	  LGCtxLeqMDI->chgainconf[0] = 1;
	  LGCtxLeqMDI->chgainconf[1] = 1;
	  LGCtxLeqMDI->chgainconf[2] = 1;
	  LGCtxLeqMDI->chgainconf[3] = 0;
	  LGCtxLeqMDI->chgainconf[4] = 1;
	  LGCtxLeqMDI->chgainconf[5] = 1;

	  LGCtxLeqMDI->chgateconf[0] = 3;
	  LGCtxLeqMDI->chgateconf[1] = 3;
	  LGCtxLeqMDI->chgateconf[2] = 3;
	  LGCtxLeqMDI->chgateconf[3] = 0;
	  LGCtxLeqMDI->chgateconf[4] = 3;
	  LGCtxLeqMDI->chgateconf[5] = 3;

	  //remember to costrain buffersizems to 400ms in case lkfs is used
	}
      else
	{
	  for (int i = 0; i < codecContext->channels; i++)
	    {
	      LGCtxLeqMDI->chgainconf[i] = 1;
	      if (i != 3)
		{
		  LGCtxLeqMDI->chgateconf[i] = 3;
		}
	      else
		{
		  LGCtxLeqMDI->chgateconf[i] = 0;
		}
	    }
	}
#endif

#ifdef SNDFILELIB
      //LGCtxLeqMDI->oldBufferBackup = malloc (sizeof (double *) * sfinfo.channels);
      LGCtxLeqMDI->oldBufferBackup = malloc (sizeof (double *) * numCPU);
      //for (int i = 0; i < sfinfo.channels; i++)
      for (int i = 0; i < numCPU; i++)
	{
	  LGCtxLeqMDI->oldBufferBackup[i] =
	    //malloc (sizeof (double) * LGCtxLeqMDI->gblocksize);
	    malloc (sizeof (double) * LGCtxLeqMDI->gblocksize *
		    sfinfo.channels);
	}
#elif defined FFMPEG
      //LGCtxLeqMDI->oldBufferBackup =
      //  malloc (sizeof (double *) * codecContext->channels);
      LGCtxLeqMDI->oldBufferBackup = malloc (sizeof (double *) * numCPU);
      //for (int i = 0; i < codecContext->channels; i++)
      for (int i = 0; i < numCPU; i++)
	{
	  *(LGCtxLeqMDI->oldBufferBackup + i) =
	    malloc (sizeof (double) * LGCtxLeqMDI->gblocksize *
		    codecContext->channels);
	}
#endif
    }				/* if (lkfs) */

  /* END LEQMDI INSERT */

#endif //this closes #ifdef DI

  if (!poly)
    {

      //End opening audio file


      eqfreqresp_db = malloc (sizeof (*eqfreqresp_db) * npoints);


      eqfreqsamples = malloc (sizeof (*eqfreqsamples) * npoints);

      eqfreqresp = malloc (sizeof (*eqfreqresp) * npoints);

      ir = malloc (sizeof (*ir) * npoints * 2);

    }				// if(!poly)
  // And what to do for floating point sample coding?
#ifdef SNDFILELIB
  switch (sfinfo.format & SF_FORMAT_SUBMASK)
    {
      // all signed bitdepth
    case 0x0001:
      bitdepth = 8;
      break;
    case 0x0002:
      bitdepth = 16;
      break;
    case 0x0003:
      bitdepth = 24;
      break;
    case 0x0004:
      bitdepth = 32;
      break;
    default:
      printf ("No known bitdepth! Exiting ...\n");
      return -1;
    }

#elif defined FFMPEG
  bitdepth = av_get_exact_bits_per_sample (codecContext->codec_id);
#ifdef DEBUG
  printf ("ffmpeg report %d bitdepth for the file.\n", bitdepth);
#endif

#endif

  if (!poly)
    {
      equalinterval3 (freqsamples, freqresp_db, eqfreqsamples, eqfreqresp_db,
		      npoints, samplingfreq, origpoints, bitdepth);
      convloglin (eqfreqresp_db, eqfreqresp, npoints);

#ifdef DEBUG
      for (int i = 0; i < npoints; i++)
	{
	  printf ("%d\t%.2f\t%.2f\t%.2f\n", i, eqfreqsamples[i],
		  eqfreqresp_db[i], eqfreqresp[i]);
	}
#endif

      inversefft2 (eqfreqresp, ir, npoints);
    }				// if (!poly)
  // read through the entire file

  struct Sum *totsum;
  totsum = malloc (sizeof (struct Sum));
  totsum->csum = 0.0;
  totsum->sum = 0.0;
  totsum->nsamples = 0;
  totsum->cmean = 0.0;
  totsum->mean = 0.0;		// Do I write anything here?
  totsum->leqm = 0.0;
  totsum->rms = 0.0;
#ifdef SNDFILELIB
  sf_count_t samples_read = 0;





#elif defined FFMPEG
  int samples_read = 0;
  int nextsample = 0;
  int dsindex = 0;
  int remaindertot = 0;
#endif

  int nchannels;
  int samplerate;
#ifdef SNDFILELIB
  nchannels = sfinfo.channels;
  samplerate = sfinfo.samplerate;
#elif defined FFMPEG
  nchannels = codecContext->channels;
  samplerate = codecContext->sample_rate;
#endif

  // Main loop through audio file

  int worker_id = 0;
  pthread_t tid[numCPU];
  struct WorkerArgs **WorkerArgsArray;
  WorkerArgsArray = malloc (sizeof (struct WorkerArgs *) * numCPU);
  int staindex = 0;		//shorttermarrayindex

#ifdef DI
  if (dolbydi)
    {
      //at present audio file must be opened two times. First time just for di to preprocess and collect dialogue information
      //do something here
      // put this later

      /* if dolbydi than buffer size will be fixed at 2048 milliseconds */
      //not necessarily it seems
      //Well my understanding is that internally do DI work like that with 75% overlapping
      size_t num_bytes;		/* Number of bytes */
      unsigned int is_eof_array[nchannels];
      unsigned int buffered_samples_array[nchannels];
      unsigned int samples_read_array[nchannels];

      /* initzialize to 0 */
      //shouldn't do it with memset?

      for (int i = 0; i < nchannels; i++)
	{
	  is_eof_array[i] = 0;
	  buffered_samples_array[i] = 0;
	  samples_read_array[i] = 0;
	}
      /* allocate the decision array */

      //allocatedidecisionarray(nchannels, numbershortperiods, shorttermdidecisionarray); //already done at this point

      /* inizialize array of DI instances, as many as channels */
      void *di_array[nchannels];
      num_bytes = di_query_mem_size (samplerate, buffersizesamplesdi);
      for (int i = 0; i < nchannels; i++)
	{
	  is_eof_array[i] = 0;
	  buffered_samples_array[i] = 0;
	  samples_read_array[i] = 0;
	  di_array[i] = (void *) malloc (num_bytes);
	}


      /* Check for errors  and configure Dialogue Intelligence */
      for (int i = 0; i < nchannels; i++)
	{
	  if (di_array[i] == NULL)
	    {
	      fprintf (stderr, "ERROR: Memory allocation failed.\n");
	      // (void) fclose(p_input);
	      //(void) fclose(p_output);
	      return ERR_CANNOT_ALLOCATE_MEMORY;
	    }
	  if (di_init
	      (di_array[i], num_bytes, samplerate, buffersizesamplesdi) < 0)
	    {
	      fprintf (stderr,
		       "ERROR: Dialogue Intelligence initialization failed.\n");
	      free (di_array[i]);
	      //(void) fclose(p_input);
	      //(void) fclose(p_output);
	      return ERR_INITIALIZATION_FAILED;

	    }
	}


      /* so many instances as channels will be running in parallel */


      //<-- START PREPROCESSING WEGEN DOLBYDI

      pthread_attr_t attr;
      pthread_attr_init (&attr);

#ifdef SNDFILELIB

      //src_data.end_of_input = 0;/* Set this later. */

      /* Start with zero to force load in while loop. */
      //src_data.input_frames = 0 ;
      //src_data.data_in = buffer ;

      //src_data.src_ratio = src_ratio ; //this could be done a separate parameter

      //src_data.data_out = src_output ;
      //src_data.output_frames = BUFFER_LEN /channels ;


      while ((samples_read =
	      sf_read_double (file, buffer, buffersizesamples)) > 0)
	{

#elif defined FFMPEG

      AVPacket readingPacket;
      av_init_packet (&readingPacket);
      int data_size = 0;
      int copiedsamples = 0;	//also pointer to position  wherein to copy into the buffer


      //int myloopcounter = 0;

      while (av_read_frame (formatContext, &readingPacket) == 0)
	{

	  //printf("External loop %d\n", myloopcounter++);


	  if (readingPacket.stream_index == audioStream->index)
	    {
	      //AVPacket decodingPacket = readingPacket;
	      int result;

	      result = avcodec_send_packet (codecContext, &readingPacket);

	      if (result < 0)
		{
		  printf ("Error submitting packet to the decoder\n");
		  exit (1);
		}

	      // Audio packets can have multiple audio frames in a single packet
	      while (readingPacket.size > 0)
		{

		  // printf("Internal 1 loop  %d\n", myloopcounter++);
		  // Try to decode the packet into a frame
		  // Some frames rely on multiple packets, so we have to make sure the frame is finished before
		  // we can use it
		  int gotFrame = 0;
		  //int result = avcodec_decode_audio4 (codecContext, frame, &gotFrame,
		  //                                    &decodingPacket);


		  result = avcodec_receive_frame (codecContext, frame);
		  if (result == 0)
		    {
		      gotFrame = 1;

		    }

		  if ((result == AVERROR (EAGAIN)) || (result == AVERROR_EOF))
		    continue;	// it was return in the ffmpeg example as it was in a function
		  else if (result < 0)
		    {
		      printf ("Error during deconding\n");
		      exit (1);
		    }



		  if (result >= 0 && gotFrame)
		    {
		      readingPacket.size -= readingPacket.size;
		      //readingPacket.data += result;
		      readingPacket.data += readingPacket.size;

		      // We now have a fully decoded audio frame
		      /*
		         #ifdef DEBUG
		         printAudioFrameInfo(codecContext, frame);
		         #endif
		       */
		      //here goes the copying to the multithreaded buffer
		      //but this buffer is in milliseconds and could be less oder more
		      //probably greater than the single frame

		      //copy as much samples as there is room in the buffer

		      if (gotFrame)
			{
			  data_size =	//this is in bytes
			    av_samples_get_buffer_size
			    (NULL,
			     codecContext->channels,
			     frame->nb_samples, codecContext->sample_fmt, 1);
			  // check if copying frame data will overflow the buffer
			  // and copy only the samples necessary to completely fill the buffer
			  // save the remaining for the next copy
			  // write a function that uses get to change data to double and copy in worker buffer

			  copiedsamples +=
			    (frame->nb_samples * frame->channels);
			  //         memcpy((char) ((void *) buffer), frame.data[0], data_size);          
			  //transfer_decoded_data(frame, WorkerArgsArray, worker_id, codecContext);
			  // nextsample is next sample of frame to be copied
			  nextsample =
			    transfer_decoded_samples (frame, buffer,
						      codecContext,
						      buffersizesamples,
						      &dsindex);
			  /*
			     #ifdef DEBUG
			     printf("Next sample index: %d\n", nextsample);
			     #endif
			   */
			  //// From here execute only if buffer is full of data

			  //if (copiedsamples >= buffersizesamples) {
			  while (copiedsamples >= buffersizesamples)
			    {
			      //printf("Internal 2 loop %d\n", myloopcounter++);
			      realnumbershortperiods++;

			      dsindex = 0;
			      //remaindertot = copiedsamples - buffersizesamples;
			      //copiedsamples = 0;

			      // store rest in another buffer
			      if (copiedsamples > buffersizesamples)
				{
				  //copiedsamples = transfer_remaining_decoded_samples(frame, remainbuffer, codecContext, nextsample, &dsindex)*frame->channels;
				  transfer_remaining_decoded_samples (frame,
								      remainbuffer,
								      codecContext,
								      nextsample,
								      &dsindex)
				    * frame->channels;
				}
			      else
				{
				  copiedsamples = 0;
				}

			      if (copiedsamples > buffersizesamples)
				{	//mistake here because I updated copiedsamples, it is no more the intended one
				  //once buffer copied, copy rest if present in buffer for next cycle
				  //but I have to update dsindex otherwise I will overwrite data on the next cycle
				  //ACHTUNG: the buffer must be the one of the next worker!
				  // ALSO: I have to account for n_samples if buffer is not full on the last buffer
				  memcpy (buffer, remainbuffer, sizeof (double) * (copiedsamples - buffersizesamples));	//here I should copy only copiedsamples-buffersizesamples
				  dsindex = copiedsamples - buffersizesamples;	//added
				  copiedsamples =
				    copiedsamples - buffersizesamples;
				  //samples_read = buffersizesamples; //pertains DI
				  //is_eof = 0; //pertains DI
				  //buffered_samples = copiedsaamples - buffersizesamples;
				  //copiedsamples already initialized to the right value
				}

			      /* Loop through channels in a short period. This loop is happening in the single threads wenn actual measuring processing takes place */
			      /* My understanding is that DI cannot use a 75% overlap of the samples, so one cannot really parallelize inside a single channel, but only between the channels. So the parallelization should be limited by the number of channels! Is it so at present? */

			      for (int ch_index = 0; ch_index < nchannels; //cannot be more than channel because per di_instance data must be coming serially
				   ch_index++)
				{
				  //printf("Internal 3 loop %d\n", myloopcounter++);
				  WorkerArgsArray[worker_id] =
				    malloc (sizeof (struct WorkerArgs));
				  WorkerArgsArray[worker_id]->
				    fullbuffer_nsamples = buffersizesamplesdi;
				  WorkerArgsArray[worker_id]->nsamples =
				    buffersizesamplesdi;
				  WorkerArgsArray[worker_id]->nch =
				    codecContext->channels;
				  WorkerArgsArray[worker_id]->npoints =
				    npoints;
				  WorkerArgsArray[worker_id]->ir = ir;
				  WorkerArgsArray[worker_id]->ptrtotsum =
				    totsum;

				  WorkerArgsArray[worker_id]->chconf =
				    channelconfcalvector;
				  /* if ((leqm10) || (dolbydi)) { */
				  WorkerArgsArray[worker_id]->shorttermindex = staindex;	// <- This was altered
				  /* if (leqm10)
				     {
				     WorkerArgsArray[worker_id]->leqm10flag = 1;  // this only if leqm10
				     WorkerArgsArray[worker_id]->shorttermarray =
				     shorttermaveragedarray;
				     //WorkerArgsArray[worker_id]->shorttermarray_di = NULL; //Do not do that just in case of leqm10 AND dolbydi
				     } */
				  /*if (dolbydi) { *///already in an if (dolbydi)
				  WorkerArgsArray[worker_id]->is_eof_array = is_eof_array;	// so later I need to reset this - only for DI
				  WorkerArgsArray[worker_id]->buffered_samples_array = buffered_samples_array;	// so later I need to reset this - only for DI
				  WorkerArgsArray[worker_id]->samples_read_array = samples_read_array;	// so later I need to reset this - only for DI
				  //WorkerArgsArray[worker_id]->leqm10flag = 0; // this only if leqm10 but do not do that just in case of multiple choises
				  //WorkerArgsArray[worker_id]->shorttermarray = NULL;  
				  WorkerArgsArray[worker_id]->
				    shorttermarray_di =
				    shorttermdidecisionarray;
				  WorkerArgsArray[worker_id]->
				    sc_shorttermarray =
				    sc_shorttermaveragedarray;
				  WorkerArgsArray[worker_id]->di_array =
				    di_array;
				  WorkerArgsArray[worker_id]->channel =
				    ch_index;
				  WorkerArgsArray[worker_id]->chgate =
				    channelgateconfvector;
				  /*} */
				  /*
				     } else {
				     WorkerArgsArray[worker_id]->shorttermindex = 0;
				     WorkerArgsArray[worker_id]->leqm10flag = 0;
				     } *///should never get here
				  /* The following will be different for preprocessing with DI and for actual measurement. Also shorter */
				  WorkerArgsArray[worker_id]->di_argbuffer =
				    malloc (sizeof (float) *
					    buffersizesamplesdi);
				  WorkerArgsArray[worker_id]->
				    samples_read_array[ch_index] =
				    (unsigned int)
				    deintsamples (WorkerArgsArray[worker_id]->
						  di_argbuffer,
						  buffersizesamplesdi, buffer, //buffersizesamplesdi, buffer,
						  ch_index, nchannels);
				  //memcpy(WorkerArgsArray[worker_id]->argbuffer, (void *)buffer, buffersizesamples*sizeof(double));
				  //pthread_attr_t attr;
				  //pthread_attr_init (&attr);
				  pthread_create (&tid[worker_id], &attr,
						  di_worker_function,
						  WorkerArgsArray[worker_id]);

				  worker_id++;


				  //   if (worker_id == numCPU) { //<-- No, see above, this cannot work like that! it must be limited by channel number not CPUs
				  if (worker_id == min(nchannels, numCPU))
				    {
				      
				      //maybe here wait for all cores to output before going on
				      // for (int idxcpu = 0; idxcpu < numCPU; idxcpu++) { //<-- Here the same!
				      for (int idxcpu = 0;
					   idxcpu < min (nchannels, numCPU);
					   idxcpu++)
					{
					  pthread_join (tid[idxcpu], NULL);
					  free (WorkerArgsArray[idxcpu]->
						di_argbuffer);
					  WorkerArgsArray[idxcpu]->
					    di_argbuffer = NULL;
					  free (WorkerArgsArray[idxcpu]);
					  WorkerArgsArray[idxcpu] = NULL;
					}
				      worker_id = 0;
				    }	//if (worker_id == nchannels)

				}	// loop through channels for preprocessing with DI

			      staindex++;	//added out of the channel loop DI

			    }	/// till here if buffer is full, but if not? It will never start doing calculations
			}	//if (gotFrame)
		      if (data_size <= 0)
			{
			  // No data yet, get more frames 
			  continue;
			}
		      // We have data, return it and come back for more later 
		      //return data_size;
		      //} //if (gotFrame)


		      /* print samples */

		      /*
		         for ( int ch_index = 0; ch_index < frame->channels; ch_index++) {
		         printf("Channel %d\n", ch_index);
		         for (int smpl_index = 0; smpl_index< frame->nb_samples; smpl_index++) {

		         printf("%0.6f\n", get(frame->data, ch_index, smpl_index, frame->channels, codecContext->sample_fmt));
		         // end print samples 
		         } //for nb_samples
		         } //for channels

		       */
		    }		//if result >= 0 && gotFrame
		  else
		    {
		      readingPacket.size = 0;
		      readingPacket.data = NULL;
		    }




		}		//while decodingPacket.size


	    }			// if readingPaket.stream

	  // You *must* call av_free_packet() after each call to av_read_frame() or else you'll leak memory
	  //av_packet_unref(&decodingPacket);
	  av_packet_unref (&readingPacket);




	  //end while worker_id
	  /// End looping cores
	}			/* while (av_read_frame(formatContext, &readingPacket) */// main loop through file



      if (worker_id != 0)
	{


	  //maybe here wait for all cores to output before going on

	  for (int idxcpu = 0; idxcpu < worker_id; idxcpu++)
	    {
	      pthread_join (tid[idxcpu], NULL);
	      free (WorkerArgsArray[idxcpu]->di_argbuffer);
	      WorkerArgsArray[idxcpu]->di_argbuffer = NULL;
	      free (WorkerArgsArray[idxcpu]);
	      WorkerArgsArray[idxcpu] = NULL;
	    }


	  worker_id = 0;




	}			//if (worker_id != 0)


      /* Adding insert for last samples at the end of the file - v. 28 */
      if (copiedsamples < buffersizesamples)
	{
	  realnumbershortperiods++;

	  for (int ch_index = 0; ch_index < nchannels; ch_index++)
	    {

	      worker_id = 0;	//this should have been already like that
	      WorkerArgsArray[worker_id] =
		malloc (sizeof (struct WorkerArgs));
	      WorkerArgsArray[worker_id]->nsamples = copiedsamples / codecContext->channels;	//is this correct?
	      WorkerArgsArray[worker_id]->fullbuffer_nsamples =
		buffersizesamplesdi;
	      WorkerArgsArray[worker_id]->nch = codecContext->channels;
	      WorkerArgsArray[worker_id]->npoints = npoints;
	      WorkerArgsArray[worker_id]->ir = ir;
	      WorkerArgsArray[worker_id]->ptrtotsum = totsum;

	      WorkerArgsArray[worker_id]->chconf = channelconfcalvector;

	      /*  if ((leqm10) || (dolbydi)) { */

	      WorkerArgsArray[worker_id]->shorttermindex = staindex;	// <- This cannot be anymore in this loop FIXME
	      /*if (leqm10)
	         {
	         WorkerArgsArray[worker_id]->leqm10flag = 1;  // this only if leqm10
	         WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
	         //WorkerArgsArray[worker_id]->shorttermarray_di = NULL;
	         } */
	      /* 
	         if (dolbydi) { *///Already in an if (dolby)
	      WorkerArgsArray[worker_id]->is_eof_array = is_eof_array;	// so later I need to reset this
	      WorkerArgsArray[worker_id]->buffered_samples_array = buffered_samples_array;	// so later I need to reset this
	      WorkerArgsArray[worker_id]->samples_read_array = samples_read_array;	// so later I need to reset this
	      //WorkerArgsArray[worker_id]->leqm10flag = 0; // this only if leqm10
	      //WorkerArgsArray[worker_id]->shorttermarray = NULL;
	      WorkerArgsArray[worker_id]->shorttermarray_di =
		shorttermdidecisionarray;
	      WorkerArgsArray[worker_id]->sc_shorttermarray =
		sc_shorttermaveragedarray;
	      WorkerArgsArray[worker_id]->di_array = di_array;
	      WorkerArgsArray[worker_id]->channel = ch_index;
	      WorkerArgsArray[worker_id]->chgate = channelgateconfvector;
	      /*
	         }
	       *//*
	         } else {
	         WorkerArgsArray[worker_id]->shorttermindex = 0;
	         WorkerArgsArray[worker_id]->leqm10flag = 0;
	         } *///but should never get here
	      dsindex = 0;
	      //remaindertot = copiedsamples - buffersizesamples;
	      //copiedsamples = 0;
	      WorkerArgsArray[worker_id]->di_argbuffer =
		malloc (sizeof (double) * buffersizesamplesdi);
	      WorkerArgsArray[worker_id]->samples_read_array[ch_index] =
		(unsigned int)
		deintsamples (WorkerArgsArray[worker_id]->di_argbuffer,
			      buffersizesamplesdi, buffer, ch_index, //buffersizesamplesdi, buffer, ch_index,
			      nchannels);
	      //   memcpy(WorkerArgsArray[worker_id]->argbuffer, (void *)buffer, copiedsamples*sizeof(double));
	      //pthread_attr_t attr;
	      //pthread_attr_init (&attr);
	      pthread_create (&tid[worker_id], &attr, di_worker_function,
			      WorkerArgsArray[worker_id]);


	      pthread_join (tid[worker_id], NULL);
	      free (WorkerArgsArray[worker_id]->di_argbuffer);
	      WorkerArgsArray[worker_id]->di_argbuffer = NULL;
	      free (WorkerArgsArray[worker_id]);
	      WorkerArgsArray[worker_id] = NULL;

	    }			// loop through channels for preprocessing with DI

	  staindex++;

	}			/*if (copiedsamples < buffersizesamples) */





      // Some codecs will cause frames to be buffered up in the decoding process. If the CODEC_CAP_DELAY flag
      // is set, there can be buffered up frames that need to be flushed, so we'll do that
      if (codecContext->codec->capabilities & AV_CODEC_CAP_DELAY)
	{
	  av_init_packet (&readingPacket);
	  // Decode all the remaining frames in the buffer, until the end is reached
	  int gotFrame = 0;
	  while (avcodec_decode_audio4
		 (codecContext, frame, &gotFrame, &readingPacket) >= 0
		 && gotFrame)
	    {
	      // We now have a fully decoded audio frame
	      // so I also need to process this
	      /*
	         #ifdef DEBUG
	         printAudioFrameInfo(codecContext, frame);
	         #endif
	       */
	    }
	  av_packet_unref (&readingPacket);
	}

#endif


#ifdef SNDFILELIB

      for (int ch_index = 0; ch_index < nchannels; ch_index++)
	{

	  WorkerArgsArray[worker_id] = malloc (sizeof (struct WorkerArgs));
	  WorkerArgsArray[worker_id]->fullbuffer_nsamples =
	    buffersizesamplesdi;
	  WorkerArgsArray[worker_id]->nsamples = buffersizesamplesdi;
	  WorkerArgsArray[worker_id]->nch = sfinfo.channels;
	  WorkerArgsArray[worker_id]->npoints = npoints;
	  WorkerArgsArray[worker_id]->ir = ir;
	  WorkerArgsArray[worker_id]->ptrtotsum = totsum;

	  WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
	  /*if (leqm10)
	     { */
	  WorkerArgsArray[worker_id]->shorttermindex = staindex;
	  /*  //FOLLOWING IF ELSE SHOULD BE ELIMINATED
	     if (leqm10) {
	     WorkerArgsArray[worker_id]->leqm10flag = 1;
	     WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
	     }      
	     else
	     {
	     //WorkerArgsArray[worker_id]->shorttermindex = 0;
	     WorkerArgsArray[worker_id]->leqm10flag = 0;
	     } */

	  WorkerArgsArray[worker_id]->is_eof_array = is_eof_array;	// so later I need to reset this - only for DI
	  WorkerArgsArray[worker_id]->buffered_samples_array = buffered_samples_array;	// so later I need to reset this - only for DI
	  WorkerArgsArray[worker_id]->samples_read_array = samples_read_array;	// so later I need to reset this - only for DI
	  //WorkerArgsArray[worker_id]->leqm10flag = 0; // this only if leqm10 but do not do that just in case of multiple choises
	  //WorkerArgsArray[worker_id]->shorttermarray = NULL;  
	  WorkerArgsArray[worker_id]->shorttermarray_di =
	    shorttermdidecisionarray;
	  WorkerArgsArray[worker_id]->sc_shorttermarray =
	    sc_shorttermaveragedarray;
	  WorkerArgsArray[worker_id]->di_array = di_array;
	  WorkerArgsArray[worker_id]->channel = ch_index;
	  WorkerArgsArray[worker_id]->chgate = channelgateconfvector;

	  WorkerArgsArray[worker_id]->di_argbuffer =	//THIS MUST BE WRONG
	    malloc (sizeof (float) * buffersizesamplesdi);
	  //   WorkerArgsArray[worder_id]->src_output = malloc(sizeof(double)*buffersizesamples); // this is for sample rate conversion, not yet used
	  //memcpy (WorkerArgsArray[worker_id]->di_argbuffer, buffer,
	  //          samples_read * sizeof (double));
	  WorkerArgsArray[worker_id]->samples_read_array[ch_index] =
	    (unsigned int)
	    deintsamples (WorkerArgsArray[worker_id]->di_argbuffer,
			  buffersizesamplesdi, buffer, ch_index, nchannels); //buffersizesamplesdi, buffer, ch_index, nchannels);
	  //pthread_attr_t attr;
	  //pthread_attr_init (&attr);
	  /*
	     if (dolbydi) { //we are already in an iif (dolbydi)
	     pthread_create(&tid[worker_id], &attr, di_worker_function, WorkerArgsArray[worker_id]);
	     } else {
	   */

	  pthread_create (&tid[worker_id], &attr, di_worker_function,
			  WorkerArgsArray[worker_id]);
	  /*
	     }
	   */
	  worker_id++;

	  //if (worker_id == numCPU) { <-- See above!
	  if (worker_id == min(nchannels,numCPU)) //this will present channel deserialization at thread level but also the possibility that numCPU < nchannels
	    {
	     
	      //maybe here wait for all cores to output before going on
	      //for (int idxcpu = 0; idxcpu < numCPU; idxcpu++) { <-- See above
	      for (int idxcpu = 0; idxcpu < min (nchannels, numCPU); idxcpu++)
		{
		  pthread_join (tid[idxcpu], NULL);
		  free (WorkerArgsArray[idxcpu]->di_argbuffer);
		  WorkerArgsArray[idxcpu]->di_argbuffer = NULL;
		  free (WorkerArgsArray[idxcpu]);
		  WorkerArgsArray[idxcpu] = NULL;
		}
	      worker_id = 0;
	    }			// if (worker_id == numCPU)


	}			// loop through channels for preprocessing with DI

      staindex++;

      //end while worker_id
      /// End looping cores
    }				// main loop through file



#endif

  //here I should wait for rest workers (< numcpu)
  //but I need to dispose of thread id.
  if (worker_id != 0)
    {				// worker_id == 0 means the number of samples was divisible through the number of cpus
      for (int idxcpu = 0; idxcpu < worker_id; idxcpu++)
	{			//worker_id is at this point one unit more than threads launched
	  pthread_join (tid[idxcpu], NULL);
	  free (WorkerArgsArray[idxcpu]->di_argbuffer);
	  WorkerArgsArray[idxcpu]->di_argbuffer = NULL;
	  free (WorkerArgsArray[idxcpu]);
	  WorkerArgsArray[idxcpu] = NULL;
	}

      worker_id = 0;

    }				// if (worker_id != 0)



  //<-- HERE ENDS PREPROCESSING THROUGH FILE

  pthread_attr_destroy (&attr);
  //free DI Instances

  for (int i = 0; i < nchannels; i++)
    {
      free (di_array[i]);
      di_array[i] = NULL;
    }

  // Some re-initializations




  staindex = 0;			//shorttermarrayindex
#ifdef SNDFILELIB
  samples_read = 0;


  //add seeking at the beginning of the file
  sf_seek (file, 0, SEEK_SET);	//never tested til now

#elif defined FFMPEG
  samples_read = 0;
  nextsample = 0;
  dsindex = 0;
  remaindertot = 0;
  //av_flush_buffer();
  // it seems this method is based on time stamps so I do not know if it will work
  av_seek_frame (formatContext, audioStream->index, 0, AVSEEK_FLAG_BACKWARD);
#endif
  if (printdiinfo)
    {
#ifdef FFMPEG
      print_di (codecContext->channels, realnumbershortperiods,
		shorttermdidecisionarray);
#elif defined SNDFILELIB
      print_di (sfinfo.channels, numbershortperiods,
		shorttermdidecisionarray);
#endif
    }
}				// if (dolbydi)

#endif //this start before "if (dolby)"

															// <--- HERE STARTS AGAIN for Leq(M)





int pthreaditer = 0;
#ifdef SNDFILELIB
double *nextworkerbufferbackup =
malloc (sizeof (double) * buffersizesamples * sfinfo.channels);
double *nextworkerbufferbackup_leqmdi =
malloc (sizeof (double) * buffersizesamples * sfinfo.channels);
#elif defined FFMPEG
double *nextworkerbufferbackup =
malloc (sizeof (double) * buffersizesamples * codecContext->channels);
double *nextworkerbufferbackup_leqmdi =
malloc (sizeof (double) * buffersizesamples * codecContext->channels);
#endif



#ifdef FFMPEG
realnumbershortperiods = 0;	// why resetting this?
#endif




pthread_attr_t attr;
pthread_attr_init (&attr);

#ifdef SNDFILELIB

															//src_data.end_of_input = 0;/* Set this later. */

															/* Start with zero to force load in while loop. */
															//     src_data.input_frames = 0 ;
															//     src_data.data_in = buffer ;

															//     src_data.src_ratio = src_ratio ; //this could be done a separate parameter

															//     src_data.data_out = src_output ;
															//     src_data.output_frames = BUFFER_LEN / sfinfo.channels ;

while ((samples_read = sf_read_double (file, buffer, buffersizesamples)) > 0)
  {



#elif defined FFMPEG


																// int myloopcounter = 0;


AVPacket readingPacket;
av_init_packet (&readingPacket);
int data_size = 0;
int copiedsamples = 0;		//also pointer to position  wherein to copy into the buffer

while (av_read_frame (formatContext, &readingPacket) == 0)
  {

    //printf("Internal 1 loop %d\n", myloopcounter++);


    if (readingPacket.stream_index == audioStream->index)
      {
	//AVPacket decodingPacket = readingPacket;
	int result;

	result = avcodec_send_packet (codecContext, &readingPacket);

	if (result < 0)
	  {
	    printf ("Error submitting packet to the decoder\n");
	    exit (1);
	  }






	// Audio packets can have multiple audio frames in a single packet
	while (readingPacket.size > 0)

	  //while (result >= 0)
	  {
	    // Try to decode the packet into a frame
	    // Some frames rely on multiple packets, so we have to make sure the frame is finished before
	    // we can use it
	    int gotFrame = 0;
	    //  int result;
	    //int result = avcodec_decode_audio4 (codecContext, frame, &gotFrame,
	    //                                        &decodingPacket);
	    //printf("Internal 2 loop %d\n", myloopcounter++);


	    result = avcodec_receive_frame (codecContext, frame);
	    if (result == 0)
	      {
		gotFrame = 1;

	      }

	    if ((result == AVERROR (EAGAIN)) || (result == AVERROR_EOF))
	      {
		printf ("Got here!!");
		continue;	//in the example would be return 
	      }
	    else if (result < 0)
	      {
		printf ("Error during deconding\n");
		exit (1);
	      }







	    if (result >= 0 && gotFrame)
	      {
		readingPacket.size -= readingPacket.size;
		//readingPacket.data += result;
		readingPacket.data += readingPacket.size;

		// We now have a fully decoded audio frame
		/*
		   #ifdef DEBUG
		   printAudioFrameInfo(codecContext, frame);
		   #endif
		 */
		//here goes the copying to the multithreaded buffer
		//but this buffer is in milliseconds and could be less oder more
		//probably greater than the single frame

		//copy as much samples as there is room in the buffer

		if (gotFrame)
		  {
		    data_size =	//this is in bytes
		      av_samples_get_buffer_size
		      (NULL,
		       codecContext->channels,
		       frame->nb_samples, codecContext->sample_fmt, 1);
		    // check if copying frame data will overflow the buffer
		    // and copy only the samples necessary to completely fill the buffer
		    // save the remaining for the next copy
		    // write a function that uses get to change data to double and copy in worker buffer

		    copiedsamples += (frame->nb_samples * frame->channels);
		    //         memcpy((char) ((void *) buffer), frame.data[0], data_size);          
		    //transfer_decoded_data(frame, WorkerArgsArray, worker_id, codecContext);
		    // nextsample is next sample of frame to be copied
		    nextsample =
		      transfer_decoded_samples (frame, buffer, codecContext,
						buffersizesamples, &dsindex);
		    /*
		       #ifdef DEBUG
		       printf("Next sample index: %d\n", nextsample);
		       #endif
		     */
		    //// From here execute only if buffer is full of data

		    //if (copiedsamples >= buffersizesamples) {
		    while (copiedsamples >= buffersizesamples)
		      {
			realnumbershortperiods++;
			WorkerArgsArray[worker_id] =
			  malloc (sizeof (struct WorkerArgs));
			WorkerArgsArray[worker_id]->nsamples =
			  buffersizesamples;
			WorkerArgsArray[worker_id]->nch =
			  codecContext->channels;
			WorkerArgsArray[worker_id]->sample_rate =
			  codecContext->sample_rate;
			WorkerArgsArray[worker_id]->npoints = npoints;
			WorkerArgsArray[worker_id]->ir = ir;
			WorkerArgsArray[worker_id]->ptrtotsum = totsum;

			WorkerArgsArray[worker_id]->chconf =
			  channelconfcalvector;
			//new
			WorkerArgsArray[worker_id]->pthread_iteration =
			  pthreaditer;
			WorkerArgsArray[worker_id]->ncpus = numCPU;
			WorkerArgsArray[worker_id]->worker_id = worker_id;
			//
			if (truepeak)
			  {
			    WorkerArgsArray[worker_id]->truepeakflag = 1;
			    WorkerArgsArray[worker_id]->truepeak =
			      truepeak_ctx;
			  }
			else
			  {
			    WorkerArgsArray[worker_id]->truepeakflag = 0;
			  }

#ifdef DI
			if ((leqm10) || (leqmlog) || (dolbydi))
			  {
#else
			if ((leqm10) || (leqmlog))
			  {
#endif
			    WorkerArgsArray[worker_id]->shorttermindex =
			      staindex;
			    WorkerArgsArray[worker_id]->leqm10flag =
			      (leqm10 ? 1 : 0);
			    WorkerArgsArray[worker_id]->leqmlogflag =
			      (leqmlog ? 1 : 0);
			    WorkerArgsArray[worker_id]->shorttermarray =
			      shorttermaveragedarray;
#ifdef DI
			    if (dolbydi)
			      {
				WorkerArgsArray[worker_id]->
				  sc_shorttermarray =
				  sc_shorttermaveragedarray;
			      }
#endif

			  }
			else
			  {
			    WorkerArgsArray[worker_id]->shorttermindex = 0;
			    WorkerArgsArray[worker_id]->leqm10flag = 0;
			    WorkerArgsArray[worker_id]->leqmlogflag = 0;

			  }

			if (lkfs)
			  {
			    WorkerArgsArray[worker_id]->Kcoeffs = coeffs;
#ifdef DI
			    if ((leqm10) || (leqmlog) || (dolbydi))
			      {
				WorkerArgsArray[worker_id]->shorttermindex =
				  staindex;
			      }
			    else
			      {
#endif
				WorkerArgsArray[worker_id]->shorttermindex =
				  staindex;
#ifdef DI
			      }
#endif
			    WorkerArgsArray[worker_id]->lg_ctx = LGCtx;
			    WorkerArgsArray[worker_id]->lkfsflag = 1;
			    WorkerArgsArray[worker_id]->lg_buffers =
			      allocateLGBuffer (LGCtx->gblocksize);
			    WorkerArgsArray[worker_id]->lg_buffers->counter =
			      staindex;
			  }
			else
			  {
			    WorkerArgsArray[worker_id]->lkfsflag = 0;
			  }

#ifdef DI
			/* LEQMDI INSERT */


			if (dolbydi)
			  {
			    WorkerArgsArray[worker_id]->shorttermindex = staindex;	// but where is this counter incremented if no lkfs or leqm10, see also preceding if
			    WorkerArgsArray[worker_id]->lg_ctx_leqmdi =
			      LGCtxLeqMDI;
			    WorkerArgsArray[worker_id]->leqmdiflag = 1;
			    WorkerArgsArray[worker_id]->lg_buffers_leqmdi =
			      allocateLGBuffer (LGCtxLeqMDI->gblocksize);
			    WorkerArgsArray[worker_id]->lg_buffers_leqmdi->
			      counter = staindex;
			  }
			else
			  {
			    WorkerArgsArray[worker_id]->leqmdiflag = 0;
			  }
			/* END LEQMDI INSERT */
#endif


			if (poly)
			  {
			    WorkerArgsArray[worker_id]->polyflag = 1;
			  }
			else
			  {
			    WorkerArgsArray[worker_id]->polyflag = 0;
			  }
			dsindex = 0;
			// store rest in another buffer
			if (copiedsamples > buffersizesamples)
			  {
			    //copiedsamples = transfer_remaining_decoded_samples(frame, remainbuffer, codecContext, nextsample, &dsindex)*frame->channels;
			    transfer_remaining_decoded_samples (frame,
								remainbuffer,
								codecContext,
								nextsample,
								&dsindex) *
			      frame->channels;
			  }
			else
			  {
			    copiedsamples = 0;
			  }
			//remaindertot = copiedsamples - buffersizesamples;
			//copiedsamples = 0;
			WorkerArgsArray[worker_id]->argbuffer =
			  malloc (sizeof (double) * buffersizesamples);
			memcpy (WorkerArgsArray[worker_id]->argbuffer,
				(void *) buffer,
				buffersizesamples * sizeof (double));
			if (lkfs)
			  {

			    memcpy (WorkerArgsArray[worker_id]->lg_ctx->
				    oldBufferBackup[worker_id],
				    (void *) nextworkerbufferbackup,
				    buffersizesamples * sizeof (double));
			    memcpy (nextworkerbufferbackup, (void *) buffer,
				    buffersizesamples * sizeof (double));
			    //            memcpy (WorkerArgsArray[worker_id]->lg_ctx->oldBufferBackup[worker_id],
			    //                  (void *) buffer, buffersizesamples * sizeof (double));

			  }
			if (dolbydi)
			  {

			    memcpy (WorkerArgsArray[worker_id]->
				    lg_ctx_leqmdi->oldBufferBackup[worker_id],
				    (void *) nextworkerbufferbackup_leqmdi,
				    buffersizesamples * sizeof (double));
			    memcpy (nextworkerbufferbackup_leqmdi, (void *) buffer,
				    buffersizesamples * sizeof (double));
			    //            memcpy (WorkerArgsArray[worker_id]->lg_ctx->oldBufferBackup[worker_id],
			    //                  (void *) buffer, buffersizesamples * sizeof (double));

			  }
			if (copiedsamples > buffersizesamples)
			  {	//mistake here because I updated copiedsamples, it is no more the intended one
			    //once buffer copied, copy rest if present in buffer for next cycle
			    //but I have to update dsindex otherwise I will overwrite data on the next cycle
			    //ACHTUNG: the buffer must be the one of the next worker!
			    // ALSO: I have to account for n_samples if buffer is not full on the last buffer
			    memcpy (buffer, remainbuffer, sizeof (double) * (copiedsamples - buffersizesamples));	//here I should copy only copiedsamples-buffersizesamples
			    dsindex = copiedsamples - buffersizesamples;	//added
			    copiedsamples = copiedsamples - buffersizesamples;
			    //copiedsamples already initialized to the right value
			  }
			//pthread_attr_t attr;
			//pthread_attr_init (&attr);
#ifdef DI
			if (dolbydi)
			  {
			    pthread_create (&tid[worker_id], &attr,
					    worker_function_gated2,
					    WorkerArgsArray[worker_id]);
			  }
			else
			  {
#endif
			    pthread_create (&tid[worker_id], &attr,
					    worker_function,
					    WorkerArgsArray[worker_id]);
#ifdef DI
			  }
#endif
			worker_id++;

/*
#ifdef DEBUG
			printf ("Worker: %d\n", worker_id);
#endif
*/


			if (worker_id == numCPU)	// add condition if buffer is full and file is at the end
			  {
			    worker_id = 0;
			    pthreaditer++;
			    //maybe here wait for all cores to output before going on
			    for (int idxcpu = 0; idxcpu < numCPU; idxcpu++)
			      {
				pthread_join (tid[idxcpu], NULL);
				free (WorkerArgsArray[idxcpu]->argbuffer);
				WorkerArgsArray[idxcpu]->argbuffer = NULL;
				if (lkfs)
				  {
				    WorkerArgsArray[idxcpu]->lg_ctx = NULL;	//here index was worker_id but it must have been wrong!
				    freeLGBuffer (WorkerArgsArray[idxcpu]->
						  lg_buffers);
				  }
				if (dolbydi)
				  {
				    WorkerArgsArray[idxcpu]->lg_ctx_leqmdi =
				      NULL;
				    freeLGBuffer (WorkerArgsArray[idxcpu]->
						  lg_buffers_leqmdi);
				  }
				free (WorkerArgsArray[idxcpu]);
				WorkerArgsArray[idxcpu] = NULL;


			      }
			  }	//if (worker_id == numCPU)

			staindex++;
		      }		/// till here if buffer is full, but if not? It will never start doing calculations
		  }		//if (gotFrame)
		if (data_size <= 0)
		  {
		    // No data yet, get more frames 
		    continue;
		  }
		// We have data, return it and come back for more later 
		//return data_size;
		//} //if (gotFrame)


		/* print samples */

		/*
		   for ( int ch_index = 0; ch_index < frame->channels; ch_index++) {
		   printf("Channel %d\n", ch_index);
		   for (int smpl_index = 0; smpl_index< frame->nb_samples; smpl_index++) {

		   printf("%0.6f\n", get(frame->data, ch_index, smpl_index, frame->channels, codecContext->sample_fmt));
		   // end print samples 
		   } //for nb_samples
		   } //for channels

		 */
	      }			//if result >= 0 && gotFrame
	    else
	      {
		readingPacket.size = 0;
		readingPacket.data = NULL;
	      }
	  }			//while decodingPacket.size

      }				// if readingPaket.stream

    // You *must* call av_free_packet() after each call to av_read_frame() or else you'll leak memory
    //av_packet_unref (&decodingPacket);
    av_packet_unref (&readingPacket);




    //end while worker_id
    /// End looping cores
  }				/* while (av_read_frame(formatContext, &readingPacket) */// main loop through file


#ifdef DEBUG

printf
  ("Situation at the end: work_id %d, copiedsamples %d, buffersizesamples %d\n",
   worker_id, copiedsamples, buffersizesamples);

#endif
																	/* Adding insert for last samples at the end of the file - v. 28 */

if (worker_id != 0)
  {



    //maybe here wait for all cores to output before going on
    for (int idxcpu = 0; idxcpu < worker_id; idxcpu++)
      {
	pthread_join (tid[idxcpu], NULL);
	free (WorkerArgsArray[idxcpu]->argbuffer);
	WorkerArgsArray[idxcpu]->argbuffer = NULL;
	if (lkfs)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx = NULL;	//here index was worker_id but it must have been wrong!
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers);
	  }
	if (dolbydi)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx_leqmdi = NULL;
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers_leqmdi);
	  }
	free (WorkerArgsArray[idxcpu]);
	WorkerArgsArray[idxcpu] = NULL;



      }

    //worker_id = 0;


  }				//if (worker_id != 0)


if (copiedsamples < buffersizesamples)
  {

    realnumbershortperiods++;
    //worker_id = 0;		//this should have been already like that
    WorkerArgsArray[worker_id] = malloc (sizeof (struct WorkerArgs));
    WorkerArgsArray[worker_id]->nsamples = copiedsamples;	//is this correct?
    WorkerArgsArray[worker_id]->nch = codecContext->channels;
    WorkerArgsArray[worker_id]->sample_rate = codecContext->sample_rate;
    WorkerArgsArray[worker_id]->npoints = npoints;
    WorkerArgsArray[worker_id]->ir = ir;
    WorkerArgsArray[worker_id]->ptrtotsum = totsum;
    //new
    WorkerArgsArray[worker_id]->pthread_iteration = pthreaditer;
    WorkerArgsArray[worker_id]->ncpus = numCPU;
    WorkerArgsArray[worker_id]->worker_id = worker_id;
    //

    WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
    if (truepeak)
      {
	WorkerArgsArray[worker_id]->truepeakflag = 1;
	WorkerArgsArray[worker_id]->truepeak = truepeak_ctx;
      }
    else
      {
	WorkerArgsArray[worker_id]->truepeakflag = 0;
      }



#ifdef DI
    if ((leqm10) || (leqmlog) || (dolbydi))
      {
#else
    if ((leqm10) || (leqmlog))
      {
#endif
	WorkerArgsArray[worker_id]->shorttermindex = staindex;
	WorkerArgsArray[worker_id]->leqm10flag = (leqm10 ? 1 : 0);
	WorkerArgsArray[worker_id]->leqmlogflag = (leqmlog ? 1 : 0);
	WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
#ifdef DI
	if (dolbydi)
	  {
	    WorkerArgsArray[worker_id]->sc_shorttermarray =
	      sc_shorttermaveragedarray;
	  }
#endif
      }
    else
      {
	WorkerArgsArray[worker_id]->shorttermindex = 0;
	WorkerArgsArray[worker_id]->leqm10flag = 0;
	WorkerArgsArray[worker_id]->leqmlogflag = 0;
      }
    if (lkfs)
      {
	WorkerArgsArray[worker_id]->Kcoeffs = coeffs;
#ifdef DI
	if ((leqm10) || (leqmlog) || (dolbydi))
	  {
	    WorkerArgsArray[worker_id]->shorttermindex = staindex;
	  }
	else
	  {
#endif

	    WorkerArgsArray[worker_id]->shorttermindex = staindex;
#ifdef DI
	  }
#endif



	WorkerArgsArray[worker_id]->lg_ctx = LGCtx;
	WorkerArgsArray[worker_id]->lkfsflag = 1;
	//WorkerArgsArray[worker_id]->lg_buffers = allocateLGBuffer(LGCtx->gblocksize);
	WorkerArgsArray[worker_id]->lg_buffers =
	  allocateLGBuffer (copiedsamples / codecContext->channels);
	WorkerArgsArray[worker_id]->lg_buffers->counter = staindex;
      }
    else
      {
	WorkerArgsArray[worker_id]->lkfsflag = 0;
      }

#ifdef DI
    /* LEQMDI INSERT */


    if (dolbydi)
      {
	WorkerArgsArray[worker_id]->shorttermindex = staindex;	// but where is this counter incremented if no lkfs or leqm10, see also preceding if
	WorkerArgsArray[worker_id]->lg_ctx_leqmdi = LGCtxLeqMDI;
	WorkerArgsArray[worker_id]->leqmdiflag = 1;
	WorkerArgsArray[worker_id]->lg_buffers_leqmdi =
	  //	  allocateLGBuffer (LGCtxLeqMDI->gblocksize); //is this causing the difference?
	  allocateLGBuffer (copiedsamples/codecContext->channels);	  
	WorkerArgsArray[worker_id]->lg_buffers_leqmdi->counter = staindex;
      }
    else
      {
	WorkerArgsArray[worker_id]->leqmdiflag = 0;
      }
    /* END LEQMDI INSERT */
#endif

    staindex++;
    if (poly)
      {
	WorkerArgsArray[worker_id]->polyflag = 1;
      }
    else
      {
	WorkerArgsArray[worker_id]->polyflag = 0;
      }
    dsindex = 0;
    //remaindertot = copiedsamples - buffersizesamples;
    //copiedsamples = 0;
    WorkerArgsArray[worker_id]->argbuffer =
      malloc (sizeof (double) * copiedsamples);
    memcpy (WorkerArgsArray[worker_id]->argbuffer, (void *) buffer,
	    copiedsamples * sizeof (double));

    if (lkfs)
      {

	memcpy (WorkerArgsArray[worker_id]->lg_ctx->
		oldBufferBackup[worker_id], (void *) nextworkerbufferbackup,
		copiedsamples * sizeof (double));
	memcpy (nextworkerbufferbackup, (void *) buffer,
		copiedsamples * sizeof (double));
	//memcpy (WorkerArgsArray[worker_id]->lg_ctx->oldBufferBackup[worker_id],
	//          (void *) buffer, copiedsamples * sizeof (double));

      }
    if (dolbydi)
      {

	memcpy (WorkerArgsArray[worker_id]->lg_ctx_leqmdi->
		oldBufferBackup[worker_id], (void *) nextworkerbufferbackup_leqmdi,
		copiedsamples * sizeof (double));
	memcpy (nextworkerbufferbackup_leqmdi, (void *) buffer,
		copiedsamples * sizeof (double));

      }



    //pthread_attr_t attr;
    //pthread_attr_init (&attr);
#ifdef DI
    if (dolbydi)
      {
	pthread_create (&tid[worker_id], &attr, worker_function_gated2,
			WorkerArgsArray[worker_id]);
      }
    else
      {
#endif

	pthread_create (&tid[worker_id], &attr, worker_function,
			WorkerArgsArray[worker_id]);
#ifdef DI
      }
#endif
    pthread_join (tid[worker_id], NULL);
    free (WorkerArgsArray[worker_id]->argbuffer);
    WorkerArgsArray[worker_id]->argbuffer = NULL;
    if (lkfs)
      {
	WorkerArgsArray[worker_id]->lg_ctx = NULL;
	freeLGBuffer (WorkerArgsArray[worker_id]->lg_buffers);
      }
    if (dolbydi)
      {
	WorkerArgsArray[worker_id]->lg_ctx_leqmdi = NULL;
	freeLGBuffer (WorkerArgsArray[worker_id]->lg_buffers_leqmdi);
      }
    free (WorkerArgsArray[worker_id]);
    WorkerArgsArray[worker_id] = NULL;



    //Store number of frames of last non full buffer
    totsum->remainder_samples = copiedsamples;


  }				/*if (copiedsamples < buffersizesamples) */





																		// Some codecs will cause frames to be buffered up in the decoding process. If the CODEC_CAP_DELAY flag
																		// is set, there can be buffered up frames that need to be flushed, so we'll do that
if (codecContext->codec->capabilities & AV_CODEC_CAP_DELAY)
  {
    av_init_packet (&readingPacket);
    // Decode all the remaining frames in the buffer, until the end is reached
    int gotFrame = 0;
    while (avcodec_decode_audio4
	   (codecContext, frame, &gotFrame, &readingPacket) >= 0 && gotFrame)
      {
	// We now have a fully decoded audio frame
	// so I also need to process this
	/*
	   #ifdef DEBUG
	   printAudioFrameInfo(codecContext, frame);
	   #endif
	 */
      }
    av_packet_unref (&readingPacket);
  }
#endif
#ifdef SNDFILELIB
WorkerArgsArray[worker_id] = malloc (sizeof (struct WorkerArgs));
WorkerArgsArray[worker_id]->nsamples = samples_read;
WorkerArgsArray[worker_id]->nch = sfinfo.channels;
WorkerArgsArray[worker_id]->sample_rate = sfinfo.samplerate;
WorkerArgsArray[worker_id]->npoints = npoints;
WorkerArgsArray[worker_id]->ir = ir;
WorkerArgsArray[worker_id]->ptrtotsum = totsum;
																		//new
WorkerArgsArray[worker_id]->pthread_iteration = pthreaditer;
WorkerArgsArray[worker_id]->ncpus = numCPU;
WorkerArgsArray[worker_id]->worker_id = worker_id;
																		//

WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
if (truepeak)
  {
    WorkerArgsArray[worker_id]->truepeakflag = 1;
    WorkerArgsArray[worker_id]->truepeak = truepeak_ctx;
  }
else
  {
    WorkerArgsArray[worker_id]->truepeakflag = 0;
  }



#ifdef DI
if ((leqm10) || (leqmlog) || (dolbydi))
  {
#else
if ((leqm10) || (leqmlog))
  {
#endif
    WorkerArgsArray[worker_id]->shorttermindex = staindex;
    WorkerArgsArray[worker_id]->leqm10flag = (leqm10 ? 1 : 0);
    WorkerArgsArray[worker_id]->leqmlogflag = (leqmlog ? 1 : 0);
    WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
#ifdef DI
    if (dolbydi)
      {
	WorkerArgsArray[worker_id]->sc_shorttermarray =
	  sc_shorttermaveragedarray;
      }
#endif
  }
else
  {
    WorkerArgsArray[worker_id]->shorttermindex = 0;
    WorkerArgsArray[worker_id]->leqm10flag = 0;
    WorkerArgsArray[worker_id]->leqmlogflag = 0;
  }

if (lkfs)
  {
    WorkerArgsArray[worker_id]->Kcoeffs = coeffs;
#ifdef DI
    if ((leqm10) || (leqmlog) || (dolbydi))
      {
#else
    if ((leqm10) || (leqmlog))
      {
#endif
	WorkerArgsArray[worker_id]->shorttermindex = staindex;
      }
    else
      {
	WorkerArgsArray[worker_id]->shorttermindex = staindex;
      }
    WorkerArgsArray[worker_id]->lg_ctx = LGCtx;
    WorkerArgsArray[worker_id]->lkfsflag = 1;
    WorkerArgsArray[worker_id]->lg_buffers =
      allocateLGBuffer (LGCtx->gblocksize);
    WorkerArgsArray[worker_id]->lg_buffers->counter = staindex;
  }
else
  {
    WorkerArgsArray[worker_id]->lkfsflag = 0;
  }


#ifdef DI
/* LEQMDI INSERT */


if (dolbydi)
  {
    WorkerArgsArray[worker_id]->shorttermindex = staindex;	// but where is this counter incremented if no lkfs or leqm10, see also preceding if
    WorkerArgsArray[worker_id]->lg_ctx_leqmdi = LGCtxLeqMDI;
    WorkerArgsArray[worker_id]->leqmdiflag = 1;
    WorkerArgsArray[worker_id]->lg_buffers_leqmdi =
      allocateLGBuffer (LGCtxLeqMDI->gblocksize);
    WorkerArgsArray[worker_id]->lg_buffers_leqmdi->counter = staindex;
  }
else
  {
    WorkerArgsArray[worker_id]->leqmdiflag = 0;
  }

																				/* END LEQMDI INSERT */
#endif


if (poly)
  {
    WorkerArgsArray[worker_id]->polyflag = 1;
  }
else
  {
    WorkerArgsArray[worker_id]->polyflag = 0;
  }

WorkerArgsArray[worker_id]->argbuffer =
malloc (sizeof (double) * buffersizesamples);
																				//   WorkerArgsArray[worder_id]->src_output = malloc(sizeof(double)*buffersizesamples); // this is for sample rate conversion, not yet used
memcpy (WorkerArgsArray[worker_id]->argbuffer, buffer,
	samples_read * sizeof (double));
if (lkfs)
  {

    memcpy (WorkerArgsArray[worker_id]->lg_ctx->oldBufferBackup[worker_id],
	    (void *) nextworkerbufferbackup, samples_read * sizeof (double));
    memcpy (nextworkerbufferbackup, (void *) buffer,
	    samples_read * sizeof (double));
    // memcpy (WorkerArgsArray[worker_id]->lg_ctx->oldBufferBackup[worker_id],
    //          buffer, samples_read * sizeof (double));

  }
if (dolbydi)
  {
    memcpy (WorkerArgsArray[worker_id]->lg_ctx_leqmdi->
	    oldBufferBackup[worker_id], (void *) nextworkerbufferbackup_leqmdi,
	    samples_read * sizeof (double));
    memcpy (nextworkerbufferbackup_leqmdi, (void *) buffer,
	    samples_read * sizeof (double));

  }


																				//pthread_attr_t attr;
																				//pthread_attr_init (&attr);
#ifdef DI
if (dolbydi)
  {
    pthread_create (&tid[worker_id], &attr, worker_function_gated2,
		    WorkerArgsArray[worker_id]);
  }
else
  {
#endif
    pthread_create (&tid[worker_id], &attr, worker_function,
		    WorkerArgsArray[worker_id]);
#ifdef DI
  }
#endif
worker_id++;


staindex++;


if (worker_id == numCPU)
  {
    worker_id = 0;
    pthreaditer++;
    //maybe here wait for all cores to output before going on
    for (int idxcpu = 0; idxcpu < numCPU; idxcpu++)
      {
	pthread_join (tid[idxcpu], NULL);
	free (WorkerArgsArray[idxcpu]->argbuffer);
	WorkerArgsArray[idxcpu]->argbuffer = NULL;
	if (lkfs)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx = NULL;
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers);
	  }
	if (dolbydi)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx_leqmdi = NULL;
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers_leqmdi);
	  }
	free (WorkerArgsArray[idxcpu]);
	WorkerArgsArray[idxcpu] = NULL;


      }
  }

																				//end while worker_id
																				// End looping cores
//Store number for frames in the last non full buffer
totsum->remainder_samples = samples_read;
}		// main loop through file //while ((samples_read = sf_read_double (file, buffer, buffersizesamples)) > 0)





																			//here I should wait for rest workers (< numcpu)
																			//but I need to dispose of thread id.
if (worker_id != 0)
  {				// worker_id == 0 means the number of samples was divisible through the number of cpus
    for (int idxcpu = 0; idxcpu < worker_id; idxcpu++)
      {				//worker_id is at this point one unit more than threads launched
	pthread_join (tid[idxcpu], NULL);
	free (WorkerArgsArray[idxcpu]->argbuffer);
	WorkerArgsArray[idxcpu]->argbuffer = NULL;
	if (lkfs)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx = NULL;
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers);
	  }
	if (dolbydi)
	  {
	    WorkerArgsArray[idxcpu]->lg_ctx_leqmdi = NULL;
	    freeLGBuffer (WorkerArgsArray[idxcpu]->lg_buffers_leqmdi);
	  }

	free (WorkerArgsArray[idxcpu]);
	WorkerArgsArray[idxcpu] = NULL;


      }

    //worker_id = 0;

  }





#endif


free (nextworkerbufferbackup);

free(nextworkerbufferbackup_leqmdi);

																			/* HERE ENDS REAL PROCESSING AFTER DI PREPROCESSING OR NOT */




																			// mean of scalar sum over duration
#ifdef DI
if (dolbydi)
  {
#ifdef FFMPEG


    dolbydifinalcomputation2 (LGCtxLeqMDI, LGCtxLeqMDI->chgateconf,
			      codecContext->channels,
			      shorttermdidecisionarray, agsthreshold);

#elif defined SNDFILELIB
    dolbydifinalcomputation2 (LGCtxLeqMDI, LGCtxLeqMDI->chgateconf,
			      sfinfo.channels,
			      shorttermdidecisionarray, agsthreshold);
#endif

  }
else if (dolbydialt)
  {

#ifdef FFMPEG
    dolbydifinalcomputation (sc_shorttermaveragedarray,
			     shorttermdidecisionarray, codecContext->channels,
			     realnumbershortperiods, tempchgate, totsum);
#elif defined SNDFILELIB

    dolbydifinalcomputation (sc_shorttermaveragedarray,
			     shorttermdidecisionarray, sfinfo.channels,
			     numbershortperiods, tempchgate, totsum);
#endif

    if (levelgated)
      {
#ifdef FFMPEG
	levelgatefinalcomputation (sc_shorttermaveragedarray,
				   leqmtosum (levelgatedthreshold, 0.000020),
				   codecContext->channels,
				   realnumbershortperiods, totsum);
#elif defined SNDFILELIB
	levelgatefinalcomputation (sc_shorttermaveragedarray,
				   leqmtosum (levelgatedthreshold, 0.000020),
				   sfinfo.channels, numbershortperiods,
				   totsum);
#endif
	printf ("Leq(M,LG): %.4f\n", totsum->lgleqm);
      }
    printf ("Leq(M,DI): %.4f\n", totsum->dgleqm);
    printf ("Program dialogue percentage is %.2f%% \n",
	    totsum->dialoguepercentual);

    meanoverduration (totsum);
    if (leqnw)
      {
	printf ("Leq(noW): %.4f\n", totsum->rms);	// Leq(no Weighting)
      }
    printf ("Leq(M): %.4f\n", totsum->leqm);
  }				// if (dolbydi) else if (dolbydialt)

#endif
meanoverduration (totsum);
if (truepeak)
  {
    printf ("True Peak Full Scale per channel:\n");
#ifdef FFMPEG
    for (int i = 0; i < codecContext->channels; i++)
      {
#elif defined SNDFILELIB
    for (int i = 0; i < sfinfo.channels; i++)
      {

#endif
	printf ("Ch %d: %.4f dBFS\n", i, log10 (truepeak_ctx->vector[i]) * 10 + 12.04);	// *10 because its power due to rectification
      }
  }
if (leqnw)
  {
    printf ("Leq(noW): %.4f\n", totsum->rms);	// Leq(no Weighting)
  }
if (lkfs)
  {
#ifdef DI
    if (dolbydi)
      {
#ifdef FFMPEG
	lkfs_finalcomputation_withdolbydi (LGCtx, LGCtx->chgateconf,
					   codecContext->channels,
					   shorttermdidecisionarray,
					   agsthreshold);
#elif defined SNDFILELIB
	lkfs_finalcomputation_withdolbydi (LGCtx, LGCtx->chgateconf,
					   sfinfo.channels,
					   shorttermdidecisionarray,
					   agsthreshold);
#endif
      }
    else
      {
#endif
#ifdef FFMPEG
	lkfs_finalcomputation (LGCtx, LGCtx->chgateconf,
			       codecContext->channels);
#elif defined SNDFILELIB
	lkfs_finalcomputation (LGCtx, LGCtx->chgateconf, sfinfo.channels);
#endif
#ifdef DI
      }
#endif
  }				// if (lkfs)
printf ("Leq(M): %.4f\n", totsum->leqm);


if (timing)
  {
    struct timespec stoptime;
    long stoptimenanoseconds;
    long executionnanoseconds;
    clock_gettime (CLOCK_MONOTONIC, &stoptime);

    if (stoptime.tv_nsec < starttime.tv_nsec)
      {
	stoptimenanoseconds = 1000000000 + stoptime.tv_nsec;
      }
    else
      {
	stoptimenanoseconds = stoptime.tv_nsec;
      }
    executionnanoseconds = stoptimenanoseconds - starttime.tv_nsec;
    printf ("Total execution time is %.6f seconds\n",
	    ((double) stoptime.tv_sec) - ((double) starttime.tv_sec) +
	    ((double) executionnanoseconds / 1000000000.00));
  }


if (leqm10)
  {
    //count shorttimeperiods through iterations!
#ifdef FFMPEG
    printf ("Number short period buffers is: %d.\n", realnumbershortperiods);
    double duration =
      ((double) realnumbershortperiods) * ((double) buffersizems / 1000.0);
#elif defined SNDFILELIB
    printf ("Number short period buffers is: %d.\n", numbershortperiods);
    double duration =
      ((double) numbershortperiods) * ((double) buffersizems / 1000.0);
#endif
    //add remainder in samples!
    //this is not precise, because the last buffer will not be full
    //Take the array with the short term accumulators
    //double interval = 10.0;
    //create a rolling average according to rolling interval
    int rollint;		// in short 10*60 = 600 sec 600/0.850 

    //how many element of the array to consider for the rollint?
    //that is how many buffersizems in the interval - interval could be parameterized(?)
    double tempint = 60.0 * longperiod / (((double) buffersizems) / 1000.0);
    rollint = (int) tempint;
    //dispose of the rest
    if (tempint - ((double) rollint) > 0)
      {
	rollint += 1;
      }


    printf ("Total duration in minutes is %.0f.\n", duration / 60.0);
    if ((duration / 60.0) < longperiod)
      {
	printf ("The audio file is too short to measure Leq(M,10m).\n");
	fclose (leqm10logfile);

	if (!leqmlog)
	  {
	    free (shorttermaveragedarray);
	    shorttermaveragedarray = NULL;
	  }

	goto skipleqm10;	//but if I really want to exit here I should free memory
      }

    //to be doubled for SNDFILELIB

    //numbershortperiods = (int) (180000.00 / (((double) codecContext->sample_rate) * (double) buffersizems/1000.00) + 1); //this is wrong, because 180000 cannot be be number of frames, also why should we devide by the sample rate if 18000 is just seconds?

    //shorttermaveragedarray = malloc(sizeof(*shorttermaveragedarray)*numbershortperiods);


#ifdef FFMPEG
    double *allenmetricarray =
      malloc (sizeof (*allenmetricarray) *
	      (realnumbershortperiods - rollint));
#elif defined SNDFILELIB
    double *allenmetricarray =
      malloc (sizeof (*allenmetricarray) * (numbershortperiods - rollint));
#endif

    //two loops
    //external loop
    int indexlong = 0;
    double temp_leqm10 = 0.0;

#ifdef FFMPEG
    while (indexlong < (realnumbershortperiods - rollint))
      {
#elif defined SNDFILELIB
    while (indexlong < (numbershortperiods - rollint))
      {
#endif
	double accumulator = 0;
	//internal loop
	double averagedaccumulator = 0;
	for (int indexshort = 0; indexshort < rollint; indexshort++)
	  {

	    accumulator += shorttermaveragedarray[indexshort + indexlong];
	  }			//end internal loop
	averagedaccumulator = accumulator / ((double) rollint);
	temp_leqm10 =
	  logleqm10 (leqm10logfile,
		     ((double) (indexlong + rollint)) *
		     ((double) buffersizems / 1000.0), averagedaccumulator);
	if (temp_leqm10 > threshold)
	  {			//See Allen article, this seems quite high...
	    allenmetricarray[indexlong] = temp_leqm10;
	  }
	else
	  {
	    allenmetricarray[indexlong] = 0.0;
	  }
	indexlong++;
      }				/* while (indexlong < < (realnumbershortperiods - rollint)) *///end external loop

    double thresholdedsum = 0;
#ifdef FFMPEG
    for (int i = 0; i < (realnumbershortperiods - rollint); i++)
      {
#elif defined SNDFILELIB
    for (int i = 0; i < (numbershortperiods - rollint); i++)
      {
#endif
	thresholdedsum += allenmetricarray[i];
      }
    //printf("Allen Metric: %d", (int) (thresholdedsum / ((double) numbershortperiods)); // But Ioan Allen seems to require minutes as unites.
    printf ("Allen metric: %d.\n", (int) (thresholdedsum / (duration / 60.0)));	// But Ioan Allen seems to require minutes as unites. But considering that the buffers are set to 750 ms it will be essentially the same, simply spreaded out times 80.
    fclose (leqm10logfile);
    if (!leqmlog)
      {
	free (shorttermaveragedarray);
	shorttermaveragedarray = NULL;
      }
    free (allenmetricarray);
    allenmetricarray = NULL;
  }



skipleqm10:

																						/*  NEW LOGLEQM */


if (leqmlog)
  {
    //count shorttimeperiods through iterations!
#ifdef FFMPEG
    printf ("Number short period buffers for Leq(M) log is: %d.\n",
	    realnumbershortperiods);
    double duration =
      ((double) realnumbershortperiods) * ((double) buffersizems / 1000.0);
#elif defined SNDFILELIB
    printf ("Number short period buffers for Leq(M) log is: %d.\n",
	    numbershortperiods);
    double duration =
      ((double) numbershortperiods) * ((double) buffersizems / 1000.0);
#endif

    printf ("Total duration of Leq(M) log  in minutes is %.0f.\n",
	    duration / 60.0);
    double accumulator = 0.0;
    double averagedaccumulator = 0.0;
    //long int samplesacc = 0;
    
#ifdef FFMPEG
    for (int i = 0; i < realnumbershortperiods; i++)
      {
#elif defined SNDFILELIB
    for (int i = 0; i < numbershortperiods; i++)
      {
#endif
	accumulator += shorttermaveragedarray[i];
	averagedaccumulator = accumulator / ((double) i + 1);
	logleqm (leqmlogfile,
		 (((double) buffersizems) * ((double) i + 1) / 1000.00),
		 averagedaccumulator);
	//samplesacc += buffersizesamples;
      }


    /* As far as all buffers are the same number of samples I should get the same result like accumulating in parallel 
       Still it would be nice to check. I should simply store the number of samples in the last buffer.
       For sndfile I can calculate this straight off at the beginning. But also for ffmpeg I should have it at the end.
     */


  /*
 * FIXME the last log in leqm log is not correct if the last buffer is not full, also there and in the check I should account for the possibility of 0 remainder
 *
 * */

#ifdef DEBUG
#ifdef FFMPEG
    double checkleqm = 10 * log10( accumulator / (((double) (realnumbershortperiods  - 1)) + ((double) totsum->remainder_samples) / buffersizesamples)) + 108.010299957;
    printf ("Buffersize in samples per channel is: %d\n", buffersizesamples / codecContext->channels);
    printf  ("Remainder samples in the last buffer are: %d\n", totsum->remainder_samples / codecContext->channels);
    printf ("Number of short period (buffers including last not full one) is: %d\n", realnumbershortperiods);
#elif defined SNDFILELIB
    double checkleqm = 10 * log10( accumulator / (((double) (numbershortperiods  - 1)) + ((double) totsum->remainder_samples) / buffersizesamples)) + 108.010299957;
    printf ("Buffersize in samples per channel is: %d\n", buffersizesamples / sfinfo.channels);
    printf  ("Remainder samples in the last buffer are: %d\n", (int) (totsum->remainder_samples / sfinfo.channels));
    printf ("Number of short period (buffers including last not full one) is: %d\n", numbershortperiods);
#endif
    printf ("Check Leq(M) from short term discrete accumulation: %.4f\n", checkleqm);

#endif

}
																							/* END NEW LOGLEQM */





#ifdef DI
if (dolbydi)
  {
    //the following is in case dolbydi was called without leqm10 
    if (shorttermaveragedarray != NULL)
      {
	free (shorttermaveragedarray);
	shorttermaveragedarray = NULL;
      }
    //allenmetricarray = NULL;
#ifdef FFMPEG
    freedidecisionarray (codecContext->channels, shorttermdidecisionarray);
    freesc_shorttermarray (codecContext->channels, sc_shorttermaveragedarray);
#elif defined SNDFILELIB
    freedidecisionarray (sfinfo.channels, shorttermdidecisionarray);
    freesc_shorttermarray (sfinfo.channels, sc_shorttermaveragedarray);

#endif
  }
#endif

if (lkfs)
  {
    /*
       #ifdef SNDFILELIB
       for (int i = 0; i < sfinfo.channels; i++)
       {
       #elif defined FFMPEG
       for (int i = 0; i < codecContext->channels; i++) */
#ifdef SNDFILELIB
    for (int i = 0; i < numCPU; i++)
      {
#elif defined FFMPEG
    for (int i = 0; i < numCPU; i++)

      {
#endif
	free (LGCtx->oldBufferBackup[i]);
	LGCtx->oldBufferBackup[i] = NULL;
      }

    free (LGCtx->oldBufferBackup);
    LGCtx->oldBufferBackup = NULL;
#ifdef SNDFILELIB
    freeLGresultarray (sfinfo.channels, LGCtx->shortperiods,
		       LGCtx->LGresultarray);
#elif defined FFMPEG
    freeLGresultarray (codecContext->channels, LGCtx->shortperiods,
		       LGCtx->LGresultarray);
#endif
    free (LGCtx->chgainconf);
    LGCtx->chgainconf = NULL;
    free (LGCtx->chgateconf);
    LGCtx->chgateconf = NULL;
    free (LGCtx);
    LGCtx = NULL;
    free (coeffs);
    coeffs = NULL;

  }

#ifdef DI
																								/* LEQMDI INSERT */

if (dolbydi)
  {
    /*
       #ifdef SNDFILELIB
       for (int i = 0; i < sfinfo.channels; i++)
       {
       #elif defined FFMPEG
       for (int i = 0; i < codecContext->channels; i++)
     */
#ifdef SNDFILELIB
    for (int i = 0; i < numCPU; i++)
      {
#elif defined FFMPEG
    for (int i = 0; i < numCPU; i++)

      {
#endif
	free (LGCtxLeqMDI->oldBufferBackup[i]);
	LGCtxLeqMDI->oldBufferBackup[i] = NULL;
      }

    free (LGCtxLeqMDI->oldBufferBackup);
    LGCtxLeqMDI->oldBufferBackup = NULL;
#ifdef SNDFILELIB
    freeLGresultarray (sfinfo.channels, LGCtxLeqMDI->shortperiods,
		       LGCtxLeqMDI->LGresultarray);
#elif defined FFMPEG
    freeLGresultarray (codecContext->channels, LGCtxLeqMDI->shortperiods,
		       LGCtxLeqMDI->LGresultarray);
#endif
    free (LGCtxLeqMDI->chgainconf);
    LGCtxLeqMDI->chgainconf = NULL;
    free (LGCtxLeqMDI->chgateconf);
    LGCtxLeqMDI->chgateconf = NULL;
    free (LGCtxLeqMDI);
    LGCtxLeqMDI = NULL;

  }


																									/* END LEQMDI INSERT */
#endif

if (truepeak)
  {
    freetruepeak (truepeak_ctx);
    truepeak_ctx = NULL;
  }



if (!poly)
  {
    free (eqfreqsamples);
    eqfreqsamples = NULL;
    free (eqfreqresp_db);
    eqfreqresp_db = NULL;
    free (eqfreqresp);
    eqfreqresp = NULL;
    free (ir);
    ir = NULL;
  }				// if (!poly)
free (channelconfcalvector);
channelconfcalvector = NULL;
free (WorkerArgsArray);
WorkerArgsArray = NULL;

free (totsum);
totsum = NULL;
free (buffer);
buffer = NULL;
#ifdef FFMPEG
free (remainbuffer);
remainbuffer = NULL;
#endif
#ifdef SNDFILELIB
sf_close (file);
#elif defined FFMPEG
av_frame_free (&frame);
avcodec_close (codecContext);
avcodec_free_context (&codecContext);
avformat_close_input (&formatContext);
#endif



pthread_mutex_destroy (&mutex);
pthread_attr_destroy (&attr);
pthread_cond_destroy (&serialsignal);
pthread_exit (NULL);

return 0;
}				// int main(...)





void *
worker_function (void *argstruct)
{

  struct WorkerArgs *thisWorkerArgs = (struct WorkerArgs *) argstruct;

  double *sumandsquarebuffer;
  double *csumandsquarebuffer;
  double *chsumaccumulator_norm;
  double *chsumaccumulator_conv;

  int copy_stepcounter;

  sumandsquarebuffer =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));


  csumandsquarebuffer =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_norm =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_conv =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));



  for (int i = 0; i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)
    {
      sumandsquarebuffer[i] = 0.0;
      csumandsquarebuffer[i] = 0.0;
      chsumaccumulator_norm[i] = 0.0;
      chsumaccumulator_conv[i] = 0.0;
    }


  for (int ch = 0; ch < thisWorkerArgs->nch; ch++)
    {

      if (thisWorkerArgs->lkfsflag)
	{


	  copy_stepcounter =
	    thisWorkerArgs->pthread_iteration * thisWorkerArgs->ncpus *
	    thisWorkerArgs->lg_ctx->subdivs +
	    thisWorkerArgs->worker_id * thisWorkerArgs->lg_ctx->subdivs;

	}

      double *normalizedbuffer;
      double *convolvedbuffer;


      normalizedbuffer =
	malloc (sizeof (double) *
		(thisWorkerArgs->nsamples / thisWorkerArgs->nch));

      convolvedbuffer =
	malloc (sizeof (double) *
		(thisWorkerArgs->nsamples / thisWorkerArgs->nch));


      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
	   n += thisWorkerArgs->nch, m++)
	{
	  // use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
	  //so not normalized but calibrated
	  normalizedbuffer[m] = thisWorkerArgs->argbuffer[n] * thisWorkerArgs->chconf[ch];	//this scale amplitude according to specified calibration

	}

      if (thisWorkerArgs->lkfsflag)
	{

	  //this is quite an interesting solution because so one avoids to read thisWorkerArgs->lg_ctx->stepcounter and consequently to devise strategies to mantain seriality of reading, which is quite difficult
	  //on the contrary it will not be necessary to have any particular precaution writing


	  for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
	       n += thisWorkerArgs->nch, m++)
	    {
	      //populate bufferA -- deinterleave
	      thisWorkerArgs->lg_buffers->bufferA[m] = thisWorkerArgs->argbuffer[n];	//
	    }




	  if (copy_stepcounter == 0)
	    {
	      //for (int n=ch, m= 0; n < thisWorkerArgs->nsamples; n += thisWorkerArgs->nch, m++) {
	      //for (int i=0; i < thisWorkerArgs->lg_ctx->gblocksize; i++) {
	      for (int i = 0;
		   i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)
		{		//this is for the case that gblocksize is incomplete like at the end of a file
		  thisWorkerArgs->lg_buffers->bufferB[i] = 0.0;	//This is for the first level gate block
		}
	    }
	  else
	    {

	      //MUTEX SHOULDN'T BE NEEDED AT THIS POINT - also copy no more needed, but deinterleaving of channels


	      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
		   n += thisWorkerArgs->nch, m++)
		{
		  thisWorkerArgs->lg_buffers->bufferB[m] =
		    thisWorkerArgs->lg_ctx->oldBufferBackup[thisWorkerArgs->
							    worker_id][n];

		}

	    }

	}			/*  <--   if (thisWorkerArgs->lkfsflag) */







      /* New True Peak */

      if (thisWorkerArgs->truepeakflag)
	{
	  double *pt_buffer;
	  if (thisWorkerArgs->lkfsflag)
	    {

	      pt_buffer = thisWorkerArgs->lg_buffers->bufferA;

	    }
	  else
	    {
	      pt_buffer = malloc (sizeof (double) * thisWorkerArgs->nsamples);

	      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
		   n += thisWorkerArgs->nch, m++)
		{
		  //populate bufferA -- deinterleave
		  pt_buffer[m] = thisWorkerArgs->argbuffer[n];	//
		}


	    }

	  pthread_mutex_lock (&mutex);
	  double temp_truepeak = thisWorkerArgs->truepeak->vector[ch];
	  pthread_mutex_unlock (&mutex);

	  temp_truepeak =
	    truepeakcheck (pt_buffer,
			   thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			   temp_truepeak,
			   thisWorkerArgs->truepeak->oversampling_ratio,
			   thisWorkerArgs->truepeak->filtertaps,
			   thisWorkerArgs->truepeak->filter_coeffs);


	  pthread_mutex_lock (&mutex);
	  thisWorkerArgs->truepeak->vector[ch] = max (temp_truepeak, thisWorkerArgs->truepeak->vector[ch]);	//THIS CORRECTION IS NOT ENOUGH that is use max() againt actual value should be enough to prevent concurrency of threads. The other possibility (not yet implemented) would be to have a vector of true peaks per thread per channel instead of a vector per channel.
	  pthread_mutex_unlock (&mutex);


	  if (!thisWorkerArgs->lkfsflag)
	    {
	      free (pt_buffer);
	      pt_buffer = NULL;
	    }
	}

      /* END INSERT */










      if (thisWorkerArgs->polyflag == 1)
	{
	  //M_filter instead of convolution
	  M_filter (convolvedbuffer, normalizedbuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch,
		    thisWorkerArgs->sample_rate);
	}
      else
	{
	  //convolution M ir
	  convolv_buff (normalizedbuffer, convolvedbuffer, thisWorkerArgs->ir,
			thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			thisWorkerArgs->npoints * 2);
	}

      //rectify, square und sum
      rectify (csumandsquarebuffer, convolvedbuffer,
	       thisWorkerArgs->nsamples / thisWorkerArgs->nch);
      rectify (sumandsquarebuffer, normalizedbuffer,
	       thisWorkerArgs->nsamples / thisWorkerArgs->nch);



      accumulatech (chsumaccumulator_norm, sumandsquarebuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch);
      accumulatech (chsumaccumulator_conv, csumandsquarebuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch);


      free (normalizedbuffer);
      normalizedbuffer = NULL;

      free (convolvedbuffer);
      convolvedbuffer = NULL;

      if (thisWorkerArgs->lkfsflag)
	{

	  int stepi = thisWorkerArgs->lg_ctx->ops;
	  //internal copy of step_counter

	  //put together samples
	  while (stepi <= thisWorkerArgs->nsamples / thisWorkerArgs->nch)
	    {			//instead of thisWorkerArgs->lg_ctx->gblocksize for block less the gblocksize at the end of file

	      /* // Where should I put this?

	         if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi) < thisWorkerArgs->lg_ctx->ops)  { // pad buffer

	         memcpy(thisWorkerArgs->lg_buffers->bufferLG, thisWorkerArgs->lg_buffers->bufferB + stepi, sizeof(double)*(thisWorkerArgs->lg_ctx->gblocksize - stepi)); //copy samples til boundary of bufferB
	         memcpy(thisWorkerArgs->lg_buffers->bufferLG + thisWorkerArgs->lg_ctx->gblocksize - stepi, thisWorkerArgs->lg_buffers->bufferA, sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi);

	         }

	       */
	      if (copy_stepcounter > 0)
		{		//that must  not be the first gate block 

		  memcpy (thisWorkerArgs->lg_buffers->bufferLG, thisWorkerArgs->lg_buffers->bufferB + stepi, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//copy samples til boundary of bufferB
		  memcpy (thisWorkerArgs->lg_buffers->bufferLG +
			  thisWorkerArgs->nsamples / thisWorkerArgs->nch -
			  stepi, thisWorkerArgs->lg_buffers->bufferA,
			  sizeof (double) * stepi);

		}
	      else
		{

		  memcpy (thisWorkerArgs->lg_buffers->bufferLG,
			  thisWorkerArgs->lg_buffers->bufferA,
			  sizeof (double) * stepi);
		  memset (thisWorkerArgs->lg_buffers->bufferLG + stepi, 0, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//Check if this correct

		  //for the case this is the first ... but preferably another solution
		}


	      /* Filtering stage one and two */
	      K_filter_stage1 (thisWorkerArgs->lg_buffers->bufferSwap,
			       thisWorkerArgs->lg_buffers->bufferLG,
			       thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			       thisWorkerArgs->Kcoeffs);
	      K_filter_stage2 (thisWorkerArgs->lg_buffers->bufferLG,
			       thisWorkerArgs->lg_buffers->bufferSwap,
			       thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			       thisWorkerArgs->Kcoeffs);
	      rectify (thisWorkerArgs->lg_buffers->bufferSwap,
		       thisWorkerArgs->lg_buffers->bufferLG,
		       thisWorkerArgs->nsamples / thisWorkerArgs->nch);
	      double tmp_sum =
		msaccumulate (thisWorkerArgs->lg_buffers->bufferSwap,
			      thisWorkerArgs->nsamples / thisWorkerArgs->nch);

	      //pthread_mutex_lock (&mutex); //THIS ALSO SHOULD NOT BE NECESSARY
	      // this should be done under mutex conditions -> shared resources!
	      // assign provisional measure of leq (not yet gated)
	      /*
	         #ifdef DEBUG
	         printf
	         ("LG block index: %d\t Channel: %d\t Sample sum: %.20f\n",
	         copy_stepcounter, ch, tmp_sum);
	         #endif
	       */
	      thisWorkerArgs->lg_ctx->LGresultarray[ch][thisWorkerArgs->shorttermindex][copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs] = tmp_sum;	//<- if below relative threshold, but postpone gating for the finale calculation
	      //pthread_mutex_unlock (&mutex);
	      copy_stepcounter++;
	      stepi += thisWorkerArgs->lg_ctx->ops;
	    }			/* while(stepi <= LGCtx->gblocksize)  */




	}			/* if (thisWorkerArgs->lkfsflag) */



    }				// loop through channels





  if (thisWorkerArgs->lkfsflag)
    {
      /* following if is to avoid incrementing stepcounter in case there where not enough samples to fill a gating block */
      if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) >=
	  thisWorkerArgs->lg_ctx->ops)
	{

	  pthread_mutex_lock (&mutex);

	  if (copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs == 0)
	    {
	      thisWorkerArgs->lg_ctx->stepcounter +=
		thisWorkerArgs->lg_ctx->subdivs;
	    }
	  else
	    {
	      thisWorkerArgs->lg_ctx->stepcounter +=
		copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs;
	    }



	  pthread_mutex_unlock (&mutex);
	}
    }
  //Create a function for this also a tag so that the worker know if he has to do this or not

  if ((thisWorkerArgs->leqm10flag) || (thisWorkerArgs->leqmlogflag))
    {
      thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex] =
	sumandshorttermavrg (chsumaccumulator_conv,
			     thisWorkerArgs->nsamples / thisWorkerArgs->nch);
#ifdef DEBUG
      printf ("%d: %.6f\n", thisWorkerArgs->shorttermindex,
	      thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex]);
#endif
    }
  pthread_mutex_lock (&mutex);
  // this should be done under mutex conditions -> shared resources!
  sumsamples (thisWorkerArgs->ptrtotsum, chsumaccumulator_norm,
	      chsumaccumulator_conv,
	      thisWorkerArgs->nsamples / thisWorkerArgs->nch);
  pthread_mutex_unlock (&mutex);


  free (sumandsquarebuffer);
  sumandsquarebuffer = NULL;

  free (csumandsquarebuffer);
  csumandsquarebuffer = NULL;

  free (chsumaccumulator_norm);
  chsumaccumulator_norm = NULL;

  free (chsumaccumulator_conv);
  chsumaccumulator_conv = NULL;

  free (thisWorkerArgs->argbuffer);
  thisWorkerArgs->argbuffer = NULL;
  // the memory pointed to by this pointer is freed in main
  // it is the same memory for all worker
  // but it is necessary to set pointer to NULL otherwise free will not work later (really?)
  thisWorkerArgs->chconf = NULL;
  pthread_exit (0);

}				//worker_function





		/* The following is to experiment with integration of Dolby DI with Leq(M) with gating like LKFS */

void *
worker_function_gated2 (void *argstruct)
{

  struct WorkerArgs *thisWorkerArgs = (struct WorkerArgs *) argstruct;

  double *sumandsquarebuffer;
  double *csumandsquarebuffer;
  double *chsumaccumulator_norm;
  double *chsumaccumulator_conv;

  int copy_stepcounter;
  int copy_stepcounter_leqmdi;

  sumandsquarebuffer =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));


  csumandsquarebuffer =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_norm =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_conv =
    malloc (sizeof (double) *
	    (thisWorkerArgs->nsamples / thisWorkerArgs->nch));



  for (int i = 0; i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)
    {
      sumandsquarebuffer[i] = 0.0;
      csumandsquarebuffer[i] = 0.0;
      chsumaccumulator_norm[i] = 0.0;
      chsumaccumulator_conv[i] = 0.0;
    }



  for (int ch = 0; ch < thisWorkerArgs->nch; ch++)
    {


      if (thisWorkerArgs->lkfsflag)
	{


	  copy_stepcounter =
	    thisWorkerArgs->pthread_iteration * thisWorkerArgs->ncpus *
	    thisWorkerArgs->lg_ctx->subdivs +
	    thisWorkerArgs->worker_id * thisWorkerArgs->lg_ctx->subdivs;

	}
      if (thisWorkerArgs->leqmdiflag)
	{


	  copy_stepcounter_leqmdi =
	    thisWorkerArgs->pthread_iteration * thisWorkerArgs->ncpus *
	    thisWorkerArgs->lg_ctx_leqmdi->subdivs +
	    thisWorkerArgs->worker_id *
	    thisWorkerArgs->lg_ctx_leqmdi->subdivs;

	}



      double *normalizedbuffer;
      double *convolvedbuffer;


      normalizedbuffer =
	malloc (sizeof (double) *
		(thisWorkerArgs->nsamples / thisWorkerArgs->nch));

      convolvedbuffer =
	malloc (sizeof (double) *
		(thisWorkerArgs->nsamples / thisWorkerArgs->nch));










      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
	   n += thisWorkerArgs->nch, m++)
	{
	  // use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
	  //so not normalized but calibrated
	  normalizedbuffer[m] = thisWorkerArgs->argbuffer[n] * thisWorkerArgs->chconf[ch];	//this scale amplitude according to specified calibration


	}


      /* lkfs */


      if (thisWorkerArgs->lkfsflag)
	{

	  //this is quite an interesting solution because so one avoids to read thisWorkerArgs->lg_ctx->stepcounter and consequently to devise strategies to mantain seriality of reading, which is quite difficult
	  //on the contrary it will not be necessary to have any particular precaution writing




	  for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
	       n += thisWorkerArgs->nch, m++)
	    {



	      //populate bufferA -- deinterleave
	      thisWorkerArgs->lg_buffers->bufferA[m] = thisWorkerArgs->argbuffer[n];	//
	    }


















	  if (copy_stepcounter == 0)
	    {			//this is problematic because lg_ctx is accessed by all threads / duplicate in lg_buffers outside the thread (before launching the thread)
	      //for (int n=ch, m= 0; n < thisWorkerArgs->nsamples; n += thisWorkerArgs->nch, m++) {
	      //for (int i=0; i < thisWorkerArgs->lg_ctx->gblocksize; i++) {
	      for (int i = 0;
		   i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)
		{		//this is for the case that gblocksize is incomplete like at the end of a file
		  thisWorkerArgs->lg_buffers->bufferB[i] = 0.0;	//This is for the first level gate block
		}
	    }
	  else
	    {
	      //MUTEX SHOULDN'T BE NEEDED AT THIS POINT - also copy no more needed, but deinterleaving of channels


	      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
		   n += thisWorkerArgs->nch, m++)
		{
		  thisWorkerArgs->lg_buffers->bufferB[m] =
		    thisWorkerArgs->lg_ctx->oldBufferBackup[thisWorkerArgs->
							    worker_id][n];
		}

	    }

	}			/*  <--   if (thisWorkerArgs->lkfsflag) */


#ifdef DI

      /* START DOLBYDILEQM  BLOCK */

      if (thisWorkerArgs->leqmdiflag)
	{

	  for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
	       n += thisWorkerArgs->nch, m++)
	    {
	      //populate bufferA -- deinterleave
	      thisWorkerArgs->lg_buffers_leqmdi->bufferA[m] = thisWorkerArgs->argbuffer[n];	//
	    }





	  if (copy_stepcounter_leqmdi == 0)
	    {			//this is problematic because lg_ctx is accessed by all threads / duplicate in lg_buffers outside the thread (before launching the thread)
	      //for (int n=ch, m= 0; n < thisWorkerArgs->nsamples; n += thisWorkerArgs->nch, m++) {
	      //for (int i=0; i < thisWorkerArgs->lg_ctx->gblocksize; i++) {
	      for (int i = 0;
		   i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)
		{		//this is for the case that gblocksize is incomplete like at the end of a file
		  thisWorkerArgs->lg_buffers_leqmdi->bufferB[i] = 0.0;	//This is for the first level gate block
		}
	    }
	  else
	    {

	      //MUTEX SHOULDN'T BE NEEDED AT THIS POINT - also copy no more needed, but deinterleaving of channels


	      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
		   n += thisWorkerArgs->nch, m++)
		{
		  thisWorkerArgs->lg_buffers_leqmdi->bufferB[m] =
		    thisWorkerArgs->lg_ctx_leqmdi->
		    oldBufferBackup[thisWorkerArgs->worker_id][n];
		}
	    }

	}			/*  <--   if (thisWorkerArgs->leqmdiflag) */

      /* DOLBYDILEQM INSERT TIL HERE */

#endif







      /* New True Peak */

      if (thisWorkerArgs->truepeakflag)
	{
	  double *pt_buffer;
	  if (thisWorkerArgs->lkfsflag)
	    {

	      pt_buffer = thisWorkerArgs->lg_buffers->bufferA;

	    }
	  else if (thisWorkerArgs->leqmdiflag)
	    {

	      pt_buffer = thisWorkerArgs->lg_buffers_leqmdi->bufferA;

	    }
	  else
	    {
	      pt_buffer = malloc (sizeof (double) * thisWorkerArgs->nsamples);

	      for (int n = ch, m = 0; n < thisWorkerArgs->nsamples;
		   n += thisWorkerArgs->nch, m++)
		{
		  //populate bufferA -- deinterleave
		  pt_buffer[m] = thisWorkerArgs->argbuffer[n];	//
		}


	    }

	  pthread_mutex_lock (&mutex);
	  double temp_truepeak = thisWorkerArgs->truepeak->vector[ch];
	  pthread_mutex_unlock (&mutex);

	  temp_truepeak =
	    truepeakcheck (pt_buffer,
			   thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			   temp_truepeak,
			   thisWorkerArgs->truepeak->oversampling_ratio,
			   thisWorkerArgs->truepeak->filtertaps,
			   thisWorkerArgs->truepeak->filter_coeffs);


	  pthread_mutex_lock (&mutex);
	  thisWorkerArgs->truepeak->vector[ch] = max (temp_truepeak, thisWorkerArgs->truepeak->vector[ch]);	//THIS CORRECTION IS NOT ENOUGH that is use max() againt actual value should be enough to prevent concurrency of threads. The other possibility (not yet implemented) would be to have a vector of true peaks per thread per channel instead of a vector per channel.
	  pthread_mutex_unlock (&mutex);


	  if (!(thisWorkerArgs->lkfsflag || thisWorkerArgs->leqmdiflag))
	    {
	      free (pt_buffer);
	      pt_buffer = NULL;
	    }
	}

      /* END INSERT */







      if (thisWorkerArgs->polyflag == 1)
	{
	  //M_filter instead of convolution
	  M_filter (convolvedbuffer, normalizedbuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch,
		    thisWorkerArgs->sample_rate);
	}
      else
	{
	  //convolution M ir
	  convolv_buff (normalizedbuffer, convolvedbuffer, thisWorkerArgs->ir,
			thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			thisWorkerArgs->npoints * 2);
	}

      //rectify, square und sum
      rectify (csumandsquarebuffer, convolvedbuffer,
	       thisWorkerArgs->nsamples / thisWorkerArgs->nch);
      rectify (sumandsquarebuffer, normalizedbuffer,
	       thisWorkerArgs->nsamples / thisWorkerArgs->nch);


      /* No more accumulating all channels, instead processing and storing results separately */


      accumulatech (chsumaccumulator_norm, sumandsquarebuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch);
      accumulatech (chsumaccumulator_conv, csumandsquarebuffer,
		    thisWorkerArgs->nsamples / thisWorkerArgs->nch);




      //Create a function for this also a tag so that the worker know if he has to do this or not

      //Why is this not present in worker_function? Because this is for each channel separately
      if (thisWorkerArgs->leqmdiflag)
	{

	  thisWorkerArgs->sc_shorttermarray[ch][thisWorkerArgs->
						shorttermindex] =
	    sumandshorttermavrg (csumandsquarebuffer,
				 thisWorkerArgs->nsamples /
				 thisWorkerArgs->nch);

	}
      /*
         #ifdef DEBUG // this must be changed to print all channels
         printf("%d: %.6f\n", thisWorkerArgs->shorttermindex, thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex]);
         #endif
       */


      free (normalizedbuffer);
      normalizedbuffer = NULL;

      free (convolvedbuffer);
      convolvedbuffer = NULL;


      /* lkfs */


      /* For LKFS according to ITU Standard:
         "The measurement interval shall be constrained such that it ends at the end of a gating block.
         Incomplete gating blocks at the end of the measurement interval are not used."

         But especially when linking against ffmpeg only here I have the exact number of frames/samples.
         So here I can have two cases: 1. nsamples / chs < stepi in this case nsamples should be discarded 
         2. stepi < nsamples / chs < gblock in  this case I would have a complete gating block and it should go through the while
       */


      /* THIS BLOCK MUST BE DOUBLED BUT ALSO SHOULD BE DEPENDENT on switch */
      if (thisWorkerArgs->lkfsflag)
	{

	  int stepi = thisWorkerArgs->lg_ctx->ops;
	  //internal copy of step_counter
	  //copy_stepcounter = thisWorkerArgs->lg_ctx->stepcounter;
	  //put together samples
	  while (stepi <= thisWorkerArgs->nsamples / thisWorkerArgs->nch)
	    {			//instead of thisWorkerArgs->lg_ctx->gblocksize for block less the gblocksize at the end of file

	      /* // Where should I put this?

	         if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi) < thisWorkerArgs->lg_ctx->ops)  { // pad buffer

	         memcpy(thisWorkerArgs->lg_buffers->bufferLG, thisWorkerArgs->lg_buffers->bufferB + stepi, sizeof(double)*(thisWorkerArgs->lg_ctx->gblocksize - stepi)); //copy samples til boundary of bufferB
	         memcpy(thisWorkerArgs->lg_buffers->bufferLG + thisWorkerArgs->lg_ctx->gblocksize - stepi, thisWorkerArgs->lg_buffers->bufferA, sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi);

	         }

	       */
	      if (copy_stepcounter > 0)
		{		//that must  not be the first gate block 

		  memcpy (thisWorkerArgs->lg_buffers->bufferLG, thisWorkerArgs->lg_buffers->bufferB + stepi, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//copy samples til boundary of bufferB
		  memcpy (thisWorkerArgs->lg_buffers->bufferLG +
			  thisWorkerArgs->nsamples / thisWorkerArgs->nch -
			  stepi, thisWorkerArgs->lg_buffers->bufferA,
			  sizeof (double) * stepi);

		}
	      else
		{

		  memcpy (thisWorkerArgs->lg_buffers->bufferLG,
			  thisWorkerArgs->lg_buffers->bufferA,
			  sizeof (double) * stepi);
		  memset (thisWorkerArgs->lg_buffers->bufferLG + stepi, 0, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//Check if this correct

		  //for the case this is the first ... but preferably another solution
		}


	      /* Filtering stage one and two */
	      K_filter_stage1 (thisWorkerArgs->lg_buffers->bufferSwap,
			       thisWorkerArgs->lg_buffers->bufferLG,
			       thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			       thisWorkerArgs->Kcoeffs);
	      K_filter_stage2 (thisWorkerArgs->lg_buffers->bufferLG,
			       thisWorkerArgs->lg_buffers->bufferSwap,
			       thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			       thisWorkerArgs->Kcoeffs);
	      rectify (thisWorkerArgs->lg_buffers->bufferSwap,
		       thisWorkerArgs->lg_buffers->bufferLG,
		       thisWorkerArgs->nsamples / thisWorkerArgs->nch);
	      double tmp_sum =
		msaccumulate (thisWorkerArgs->lg_buffers->bufferSwap,
			      thisWorkerArgs->nsamples / thisWorkerArgs->nch);

	      //pthread_mutex_lock (&mutex); //NOT NECESSARY ANYMORE
	      // this should be done under mutex conditions -> shared resources!
	      // assign provisional measure of leq (not yet gated)
	     /* 
	         #ifdef DEBUG
	         printf
	         ("LG block index: %d\t Channel: %d\t Sample sum: %.20f\n",
	         copy_stepcounter, ch, tmp_sum);
	         #endif
	       */
	      thisWorkerArgs->lg_ctx->LGresultarray[ch][thisWorkerArgs->shorttermindex][copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs] = tmp_sum;	//<- if below relative threshold, but postpone gating for the finale calculation
	      //pthread_mutex_unlock (&mutex);
	      copy_stepcounter++;
	      stepi += thisWorkerArgs->lg_ctx->ops;
	    }			/* while(stepi <= LGCtx->gblocksize)  */


	}			/* if (thisWorkerArgs->lkfsflag) */


#ifdef DI

      /* LEQMDI INSERT 2 START HERE  */
      if (thisWorkerArgs->leqmdiflag)
	{

	  int stepi = thisWorkerArgs->lg_ctx_leqmdi->ops;
	  //internal copy of step_counter
	  //copy_stepcounter = thisWorkerArgs->lg_ctx_leqmdi->stepcounter;
	  //put together samples
	  while (stepi <= thisWorkerArgs->nsamples / thisWorkerArgs->nch)
	    {			//instead of thisWorkerArgs->lg_ctx->gblocksize for block less the gblocksize at the end of file

	      /* // Where should I put this?

	         if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi) < thisWorkerArgs->lg_ctx->ops)  { // pad buffer

	         memcpy(thisWorkerArgs->lg_buffers->bufferLG, thisWorkerArgs->lg_buffers->bufferB + stepi, sizeof(double)*(thisWorkerArgs->lg_ctx->gblocksize - stepi)); //copy samples til boundary of bufferB
	         memcpy(thisWorkerArgs->lg_buffers->bufferLG + thisWorkerArgs->lg_ctx->gblocksize - stepi, thisWorkerArgs->lg_buffers->bufferA, sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch) - stepi);

	         }

	       */
	      if (copy_stepcounter_leqmdi > 0)
		{		//that must  not be the first gate block 

		  memcpy (thisWorkerArgs->lg_buffers_leqmdi->bufferLG, thisWorkerArgs->lg_buffers_leqmdi->bufferB + stepi, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//copy samples til boundary of bufferB
		  memcpy (thisWorkerArgs->lg_buffers_leqmdi->bufferLG +
			  thisWorkerArgs->nsamples / thisWorkerArgs->nch -
			  stepi, thisWorkerArgs->lg_buffers_leqmdi->bufferA,
			  sizeof (double) * stepi);

		}
	      else
		{

		  memcpy (thisWorkerArgs->lg_buffers_leqmdi->bufferLG,
			  thisWorkerArgs->lg_buffers_leqmdi->bufferA,
			  sizeof (double) * stepi);
		  memset (thisWorkerArgs->lg_buffers_leqmdi->bufferLG + stepi, 0, sizeof (double) * (thisWorkerArgs->nsamples / thisWorkerArgs->nch - stepi));	//Check if this correct

		  //for the case this is the first ... but preferably another solution
		}

	      /*

	         - SCALING ACCORDING TO M FILTER MUST BE ADDED

	       */


	      for (int i = 0;
		   i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++)

		{
		  // use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
		  //so not normalized but calibrated
		  thisWorkerArgs->lg_buffers_leqmdi->bufferSwap[i] = thisWorkerArgs->lg_buffers_leqmdi->bufferLG[i] * thisWorkerArgs->chconf[ch];	//this scale amplitude according to specified calibration

		}






	      if (thisWorkerArgs->polyflag == 1)
		{
		  //M_filter instead of convolution
		  M_filter (thisWorkerArgs->lg_buffers_leqmdi->bufferLG,
			    thisWorkerArgs->lg_buffers_leqmdi->bufferSwap,
			    thisWorkerArgs->nsamples / thisWorkerArgs->nch,
			    thisWorkerArgs->sample_rate);
		}
	      else
		{
		  //convolution M ir
		  convolv_buff (thisWorkerArgs->lg_buffers_leqmdi->bufferSwap,
				thisWorkerArgs->lg_buffers_leqmdi->bufferLG,
				thisWorkerArgs->ir,
				thisWorkerArgs->nsamples /
				thisWorkerArgs->nch,
				thisWorkerArgs->npoints * 2);
		}

	      rectify (thisWorkerArgs->lg_buffers_leqmdi->bufferSwap,
		       thisWorkerArgs->lg_buffers_leqmdi->bufferLG,
		       thisWorkerArgs->nsamples / thisWorkerArgs->nch);
	      double tmp_sum =
		msaccumulate (thisWorkerArgs->lg_buffers_leqmdi->bufferSwap,
			      thisWorkerArgs->nsamples / thisWorkerArgs->nch);
	      //pthread_mutex_lock (&mutex); //NO MORE NECESSARY
	      // this should be done under mutex conditions -> shared resources!
	      // assign provisional measure of leq (not yet gated)
	     /* 
	         #ifdef DEBUG
	         printf
	         ("LG Leq(M,DI) block index: %d\t Channel: %d\t Sample sum: %.20f\n",
	         copy_stepcounter_leqmdi, ch, tmp_sum);
	         #endif
	       */
	      thisWorkerArgs->lg_ctx_leqmdi->LGresultarray[ch][thisWorkerArgs->shorttermindex][copy_stepcounter_leqmdi % thisWorkerArgs->lg_ctx_leqmdi->subdivs] = tmp_sum;	//<- if below relative threshold, but postpone gating for the finale calculation
	      //pthread_mutex_unlock (&mutex);
	      copy_stepcounter_leqmdi++;
	      stepi += thisWorkerArgs->lg_ctx_leqmdi->ops;
	    }			/* while(stepi <= LGCtx->gblocksize)  */


	}			/* if (thisWorkerArgs->leqmdiflag) */

      /* END LEQMDI INSERT */
#endif

    }				// loop through channels



  if (thisWorkerArgs->lkfsflag)
    {
      /* following if is to avoid incrementing stepcounter in case there where not enough samples to fill a gating block */
      if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) >=
	  thisWorkerArgs->lg_ctx->ops)
	{
	  pthread_mutex_lock (&mutex);
	  if (copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs == 0)
	    {
	      thisWorkerArgs->lg_ctx->stepcounter +=
		thisWorkerArgs->lg_ctx->subdivs;
	    }
	  else
	    {
	      thisWorkerArgs->lg_ctx->stepcounter +=
		copy_stepcounter % thisWorkerArgs->lg_ctx->subdivs;
	    }

	  pthread_mutex_unlock (&mutex);
	}
    }


  /* LEQMDI INSERT 3 */


#ifdef DI

  if (thisWorkerArgs->leqmdiflag)
    {
      /* following if is to avoid incrementing stepcounter in case there where not enough samples to fill a gating block */
      if ((thisWorkerArgs->nsamples / thisWorkerArgs->nch) >=
	  thisWorkerArgs->lg_ctx_leqmdi->ops)
	{
	  pthread_mutex_lock (&mutex);
	  if (copy_stepcounter_leqmdi %
	      thisWorkerArgs->lg_ctx_leqmdi->subdivs == 0)
	    {
	      thisWorkerArgs->lg_ctx_leqmdi->stepcounter +=
		thisWorkerArgs->lg_ctx_leqmdi->subdivs;
	    }
	  else
	    {
	      thisWorkerArgs->lg_ctx_leqmdi->stepcounter +=
		copy_stepcounter_leqmdi %
		thisWorkerArgs->lg_ctx_leqmdi->subdivs;
	    }

	  pthread_mutex_unlock (&mutex);
	}
    }


  /* END LEQMDI INSERT 3 */
#endif






  if ((thisWorkerArgs->leqm10flag) || (thisWorkerArgs->leqmlogflag))
    {
      thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex] =
	sumandshorttermavrg (chsumaccumulator_conv,
			     thisWorkerArgs->nsamples / thisWorkerArgs->nch);
#ifdef DEBUG
      printf ("%d: %.6f\n", thisWorkerArgs->shorttermindex,
	      thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex]);
#endif
    }








  pthread_mutex_lock (&mutex);
  // this should be done under mutex conditions -> shared resources!
  sumsamples (thisWorkerArgs->ptrtotsum, chsumaccumulator_norm,
	      chsumaccumulator_conv,
	      thisWorkerArgs->nsamples / thisWorkerArgs->nch);
  pthread_mutex_unlock (&mutex);


  free (sumandsquarebuffer);
  sumandsquarebuffer = NULL;

  free (csumandsquarebuffer);
  csumandsquarebuffer = NULL;

  free (chsumaccumulator_norm);
  chsumaccumulator_norm = NULL;

  free (chsumaccumulator_conv);
  chsumaccumulator_conv = NULL;

  free (thisWorkerArgs->argbuffer);
  thisWorkerArgs->argbuffer = NULL;
  // the memory pointed to by this pointer is freed in main
  // it is the same memory for all worker
  // but it is necessary to set pointer to NULL otherwise free will not work later (really?)
  thisWorkerArgs->chconf = NULL;
  pthread_exit (0);

}				//worker_function_gated2







#ifdef FFMPEG
int
transfer_decoded_data (AVFrame * ptr_frame,
		       struct WorkerArgs **dptrWorkerArgs, int w_id,
		       AVCodecContext * codecCon)
{
  int doublesample_index = 0;	//this is to index the millisecond buffer of the worker
  for (int ch_index = 0; ch_index < ptr_frame->channels; ch_index++)
    {

      for (int smpl_index = 0; smpl_index < ptr_frame->nb_samples;
	   smpl_index++)
	{			//limit to buffer measure and return rest?

	  dptrWorkerArgs[w_id]->argbuffer[doublesample_index++] =
	    get (ptr_frame->data, ch_index, smpl_index, ptr_frame->channels,
		 codecCon->sample_fmt);
	  // end print samples 
	}			//for nb_samples
    }				//for channels
  return ptr_frame->nb_samples;
}

						// will return 0 or the index of next to be copied sample, if not enough room in buffer
int
transfer_decoded_samples (AVFrame * ptr_frame, double *buf,
			  AVCodecContext * codecCon, int buffersizs,
			  int *doublesample_index)
{
  // static int doublesample_index = 0; //this is to index the millisecond buffer of the worker
  //Yes the sample loop must be the external one, as I work with interleaved channels, not planar channels
  for (int smpl_index = 0; smpl_index < ptr_frame->nb_samples; smpl_index++)
    {				//limit to buffer measure and return rest?
      for (int ch_index = 0; ch_index < ptr_frame->channels; ch_index++)
	{
	  buf[(*doublesample_index)++] =
	    get (ptr_frame->data, ch_index, smpl_index, ptr_frame->channels,
		 codecCon->sample_fmt);
	  /*
	     #ifdef DEBUG
	     printf("%0.5f\n", buf[(*doublesample_index)-1]);
	     #endif
	   */
	  if ((*doublesample_index) == buffersizs)
	    {
	      (*doublesample_index) = 0;
	      //return (ptr_frame->nb_samples - smpl_index) * ptr_frame->channels;
	      return smpl_index + 1;
	    }
	  // end print samples 
	}			//for nb_samples
    }				//for channels
  return ptr_frame->nb_samples;
}


int
transfer_remaining_decoded_samples (AVFrame * ptr_frame, double *bufremain,
				    AVCodecContext * codecCon, int nxtsmpl,
				    int *doublesample_index)
{
  //int doublesample_index = 0; //this is to index the millisecond buffer of the worker
  //Yes the sample loop must be the external one, as I work with interleaved channels, not planar channels
  for (int smpl_index = nxtsmpl; smpl_index < ptr_frame->nb_samples;
       smpl_index++)
    {				//limit to buffer measure and return rest?
      for (int ch_index = 0; ch_index < ptr_frame->channels; ch_index++)
	{
	  bufremain[(*doublesample_index)++] =
	    get (ptr_frame->data, ch_index, smpl_index, ptr_frame->channels,
		 codecCon->sample_fmt);
	  /*
	     #ifdef DEBUG
	     printf("%0.5f\n", bufremain[(*doublesample_index)-1]);
	     #endif                     
	   */
	  // end print samples 
	}			//for nb_samples
    }				//for channels
  return (ptr_frame->nb_samples) - nxtsmpl;	//but what if also the second round at the same frame fill the buffer?
}

#endif

						//to get impulse response frequency response at equally spaced intervals is needed

int
equalinterval (double *freqsamples, double *freqresp, double *eqfreqsamples,
	       double *eqfreqresp, int points, int samplingfreq,
	       int origpoints)
{
  double freq;
  // int findex = 0;
  // int rindex = 0;
  double pass = ((double) (samplingfreq >> 1)) / ((double) points);
  for (int ieq = 0, i = 0; ieq < points; ieq++)
    {
      freq = ieq * pass;
      eqfreqsamples[ieq] = freq;

      if ((freq == 0.0) || (freq < freqsamples[1]))
	{
	  eqfreqresp[ieq] = freqresp[0];
	  continue;
	}
      else
	{

	  if ((freq >= freqsamples[i]) && (freq < freqsamples[i + 1]))
	    {
	      eqfreqresp[ieq] =
		((freqresp[i + 1] - freqresp[i]) / (freqsamples[i + 1] -
						    freqsamples[i])) * (freq -
									freqsamples
									[i]) +
		freqresp[i];
	    }
	  else if (freq >= freqsamples[i + 1])
	    {
	      while (freq >= freqsamples[i + 1])
		{
		  i++;
		  if ((i + 1) >= origpoints)
		    {
		      break;
		    }
		}
	      if ((i + 1) < origpoints)
		{
		  eqfreqresp[ieq] =
		    ((freqresp[i + 1] - freqresp[i]) / (freqsamples[i + 1] -
							freqsamples[i])) *
		    (freq - freqsamples[i]) + freqresp[i];
		}
	      else
		{
		  eqfreqresp[ieq] =
		    ((1 - freqresp[i]) / (((double) (samplingfreq >> 1)) -
					  freqsamples[i])) * (freq -
							      freqsamples[i])
		    + freqresp[i];
		}
	    }
	}
    }
  return 0;
}





						//the following is different from version 1 because interpolate between db and not linear. Conversion from db to lin must be done after.
						//it is also different for the way it interpolates between DC and 31 Hz
						// Pay attention that also arguments to the functions are changed
int
equalinterval2 (double freqsamples[], double freqresp_db[],
		double *eqfreqsamples, double *eqfreqresp, int points,
		int samplingfreq, int origpoints, int bitdepthsoundfile)
{
  double freq;
  //calculate miminum attenuation depending on the bitdeph (minus one), that is −6.020599913 dB per bit in eccess to sign
  double dcatt = ((double) (bitdepthsoundfile - 1)) * (-6.020599913);	// + 20.00; //in dB
  //double dcatt = -90.3;
  double pass = ((double) (samplingfreq >> 1)) / ((double) points);
  for (int ieq = 0, i = 0; ieq < points; ieq++)
    {
      freq = ieq * pass;
      eqfreqsamples[ieq] = freq;
      if (freq == 0.0)
	{			// I would like to cange this
	  eqfreqresp[ieq] = dcatt;
	}
      else if (freq < freqsamples[0])
	{			// this has a lot of influence on final Leq(M) value
	  eqfreqresp[ieq] =
	    ((freqresp_db[0] - dcatt) / (freqsamples[0] - 0)) * freq + dcatt;
	  //eqfreqresp[ieq] = freqresp_db[0]; // Is this meaningful? Shouldn't I interpolate between 0 Hz and 31 Hz? Otherwise for DC I have -35.5 dB
	  continue;
	}
      else
	{
	  if ((freq >= freqsamples[i]) && (freq < freqsamples[i + 1]))
	    {
	      eqfreqresp[ieq] =
		((freqresp_db[i + 1] - freqresp_db[i]) / (freqsamples[i + 1] -
							  freqsamples[i])) *
		(freq - freqsamples[i]) + freqresp_db[i];
	    }
	  else if (freq >= freqsamples[i + 1])
	    {
	      while (freq >= freqsamples[i + 1])
		{
		  i++;
		  if ((i + 1) >= origpoints)
		    {
		      break;
		    }
		}
	      if ((i + 1) < origpoints)
		{
		  eqfreqresp[ieq] =
		    ((freqresp_db[i + 1] -
		      freqresp_db[i]) / (freqsamples[i + 1] -
					 freqsamples[i])) * (freq -
							     freqsamples[i]) +
		    freqresp_db[i];
		}
	      else
		{
		  eqfreqresp[ieq] =
		    ((1 - freqresp_db[i]) / (((double) (samplingfreq >> 1)) -
					     freqsamples[i])) * (freq -
								 freqsamples
								 [i]) +
		    freqresp_db[i];
		}
	    }
	}
    }
  return 0;
}

int
equalinterval3 (double freqsamples[], double freqresp_db[],
		double *eqfreqsamples, double *eqfreqresp, int points,
		int samplingfreq, int origpoints, int bitdepthsoundfile)
{
  double freq;
  double firstslope =
    (freqresp_db[1] - freqresp_db[0]) / (freqsamples[1] - freqsamples[0]);
  /* DC attenuation is calculated from the first slope of attenuation curve sampling points */
  double dcatt = freqresp_db[0] - firstslope * freqsamples[0];
  double pass = ((double) (samplingfreq >> 1)) / ((double) points);
  for (int ieq = 0, i = 0; ieq < points; ieq++)
    {
      freq = ieq * pass;
      eqfreqsamples[ieq] = freq;
      if (freq < freqsamples[0])
	{			/* If frequency is less than first sampling point of M definition, slope of first and second point is used */// this has a lot of influence on final Leq(M) value
	  eqfreqresp[ieq] = firstslope * freq + dcatt;
	  //eqfreqresp[ieq] = freqresp_db[0]; // Is this meaningful? Shouldn't I interpolate between 0 Hz and 31 Hz? Otherwise for DC I have -35.5 dB
	  continue;
	}
      else
	{			/* all frequencies in between sampling frequencies of M definition */
	  if ((freq >= freqsamples[i]) && (i == (origpoints - 1)))
	    {
	      eqfreqresp[ieq] =
		((freqresp_db[origpoints - 1] -
		  freqresp_db[origpoints - 2]) / (freqsamples[origpoints -
							      1] -
						  freqsamples[origpoints -
							      2])) * (freq -
								      freqsamples
								      [origpoints
								       - 1]) +
		freqresp_db[origpoints - 1];
	    }
	  else if ((freq >= freqsamples[i]) && (freq < freqsamples[i + 1]))
	    {
	      eqfreqresp[ieq] =
		((freqresp_db[i + 1] - freqresp_db[i]) / (freqsamples[i + 1] -
							  freqsamples[i])) *
		(freq - freqsamples[i]) + freqresp_db[i];
	    }
	  else if (freq >= freqsamples[i + 1])
	    {
	      while (freq >= freqsamples[i + 1])
		{
		  i++;
		  if ((i + 1) >= origpoints)
		    {
		      break;
		    }
		}		// while
	      if ((i + 1) < origpoints)
		{
		  eqfreqresp[ieq] =
		    ((freqresp_db[i + 1] -
		      freqresp_db[i]) / (freqsamples[i + 1] -
					 freqsamples[i])) * (freq -
							     freqsamples[i]) +
		    freqresp_db[i];
		}
	      else
		{
		  eqfreqresp[ieq] =
		    ((1 - freqresp_db[i]) / (((double) (samplingfreq >> 1)) -
					     freqsamples[i])) * (freq -
								 freqsamples
								 [i]) +
		    freqresp_db[i];
		}
	    }
	}
    }
  return 0;
}





int
convloglin (double *in, double *out, int points)
{
  for (int i = 0; i < points; i++)
    {
      out[i] = pow (10, (in[i] / 20.0));
    }

  return 0;
}


double
convlinlog_single (double in)
{
  double out;
  out = log10 (in) * 20.0f;
  return out;
}


double
convloglin_single (double in)
{
  double out;
  out = powf (10, in / 20.0f);
  return out;
}

						// convolution

int
convolv_buff (double *sigin, double *sigout, double *impresp, int sigin_dim,
	      int impresp_dim)
{


  double sum = 0.0;
  for (int i = 0; i < sigin_dim; i++)
    {

      int m = i;
      for (int l = impresp_dim - 1; l >= 0; l--, m++)
	{
	  if (m >= sigin_dim)
	    {
	      m -= sigin_dim;
	    }
	  sum += sigin[m] * impresp[l];
	}
      sigout[i] = sum;
      sum = 0.0;
    }
  return 0;

}


void
inversefft2 (double *eqfreqresp, double *ir, int npoints)
{
  for (int n = 0; n < npoints; n++)
    {
      double parsum = 0.0;
      double partial = 0.0;

      for (int m = 1; m <= npoints - 1; m++)
	{
	  partial =
	    cos (2.0 * M_PI * ((double) m) *
		 ((((double) n) -
		   (((double) npoints) * 2.0 -
		    1) / 2) / (((double) npoints) * 2.0)));
	  parsum = parsum + eqfreqresp[m] * partial;
	}
      ir[n] = (eqfreqresp[0] + 2.0 * parsum) / ((double) npoints * 2.0);
      /*
         #ifdef DEBUG
         printf("%.4f\n", ir[n]);
         #endif
       */
    }
  for (int n = 0; n < npoints; n++)
    {
      ir[npoints + n] = ir[npoints - (n + 1)];
      /*
         #ifdef DEBUG
         printf("%.4f\n", ir[npoints+n]);
         #endif
       */
    }


}

						// scale input according to required calibration
						// this could be different for certain digital cinema formats
double
inputcalib (double dbdiffch)
{

  double coeff = pow (10, dbdiffch / 20);
  return coeff;

}

						//Change this in a Macro for speed 
						//rectify, square and sum
int
rectify (double *squared, double *inputsamples, int nsamples)
{
  for (int i = 0; i < nsamples; i++)
    {
      squared[i] = (double) powf (inputsamples[i], 2);
    }
  return 0;

}

int
initbuffer (double *buffertoinit, int nsamples)
{
  for (int i = 0; i < nsamples; i++)
    {
      buffertoinit[i] = 0.0;

    }
  return 0;
}

int
accumulatech (double *chaccumulator, double *inputchannel, int nsamples)
{
  for (int i = 0; i < nsamples; i++)
    {
      chaccumulator[i] += inputchannel[i];
    }
  return 0;
}

double
msaccumulate (double *inputbuffer, int nsamples)
{
  double sum = 0.0;
  for (int i = 0; i < nsamples; i++)
    {
      sum += inputbuffer[i];
    }
  return (sum / ((double) nsamples));
}

#ifdef DI

int
accumulatechwithgate (double *chaccumulator, double *inputchannel,
		      int nsamples, int chgateconf,
		      struct WorkerArgs *workerargsinstance,
		      int channel_index)
{
  /* chgate == 0 -> no gate
     chgate == 1 -> level gate
     chgate == 2 -> dialogue gate
   */
  if (chgateconf == 2)
    {				//if channel is configured for Dolby DI 
      for (int i = 0; i < nsamples; i++)
	{
	  if ((workerargsinstance->shorttermarray_di)[channel_index]
	      [workerargsinstance->shorttermindex] == 1)
	    {
	      chaccumulator[i] += inputchannel[i];
	    }
	}
    }
  return 0;
}

#endif


int
sumsamples (struct Sum *ts, double *inputsamples, double *cinputsamples,
	    int nsamples)
{
  ts->nsamples += nsamples;
  for (int i = 0; i < nsamples; i++)
    {
      ts->sum += inputsamples[i];
      ts->csum += cinputsamples[i];
    }
  return 0;

}

int
meanoverduration (struct Sum *oldsum)
{
  oldsum->mean = pow (oldsum->sum / ((double) oldsum->nsamples), 0.500);
  oldsum->cmean = pow (oldsum->csum / ((double) oldsum->nsamples), 0.500);
  oldsum->rms = 20 * log10 (oldsum->mean) + 108.010299957;
  if (oldsum->rms < 0.0)
    {
      oldsum->rms = 0.0;
    }
  oldsum->leqm = 20 * log10 (oldsum->cmean) + 108.010299957;	//
  if (oldsum->leqm < 0.0)
    {
      oldsum->leqm = 0.0;
    }

  /*
     How the final offset is calculated without reference to a test tone:
     P0 is the SPL reference 20 uPa

     Reference SPL is RMS ! So 85 SPL over 20 uPa is 10^4.25 x 0.000020 = 0.355655882 Pa (RMS), 
     but Peak value is 0.355655882 x sqr(2) = 0.502973372 that is 20 x log ( 0.502973372 / 0.000020) = 88.010299957
     To that one has to add the 20 dB offset of the reference -20dBFS: 88.010299957 + 20.00 = 108.010299957 

     Well, this interpretation/explanation is somehow a misconception, because be it that you use a RMS, an average or a peak meter,
     at the end you should use the same meter for both reference as well as the lower threshold 20uPa (that is even 20uPa
     should be measured with the same sort of meter). This will result in a reading offset of ca. -3 dB. This is also what
     I had implemented.
   */
  /*But ISO 21727:2004(E) ask for a reference level "measured using an average responding meter". So reference level is not 0.707, but 0.637 = 2/pi
   */
  return 0;
}


double
sumandshorttermavrg (double *channelaccumulator, int nsamples)
{
  double stsum = 0.0;
  for (int i = 0; i < nsamples; i++)
    {
      stsum += channelaccumulator[i];

    }
  return stsum / (double) nsamples;
}

void
logleqm (FILE * filehandle, double featuretimesec, double temp_leqm)
{
  temp_leqm = 20 * log10 (pow (temp_leqm, 0.500)) + 108.010299957;
  if (temp_leqm < 0.0)
    {
      temp_leqm = 0.0;
    }
  fprintf (filehandle, "%.4f", featuretimesec);
  fprintf (filehandle, "\t");
  fprintf (filehandle, "%.4f\n", temp_leqm);


}

double
logleqm10 (FILE * filehandle, double featuretimesec, double longaverage)
{
  double leqm10 = 20 * log10 (pow (longaverage, 0.500)) + 108.010299957;
  if (leqm10 < 0.0)
    {
      leqm10 = 0.0;
    }
  fprintf (filehandle, "%.4f", featuretimesec);
  fprintf (filehandle, "\t");
  fprintf (filehandle, "%.4f\n", leqm10);
  return leqm10;
}

#ifdef DI

void
savedidecision (uint8_t ** dibytearray, int shortperiodidx, uint8_t dibyte,
		int chnumb)
{
  //dibytearray[chnumb][shortperiodidx] |= (dibyte << chnumb);
  dibytearray[chnumb][shortperiodidx] = dibyte;
}

#endif

double ***
allocateLGresultarray (int nchannels, int shortperiods, int olsteps)
{
  //i channels, j overlap step, k shortperiods
  double ***pt_LGresults = malloc (sizeof (double **) * nchannels);
  for (int i = 0; i < nchannels; i++)
    {
      *(pt_LGresults + i) = malloc (sizeof (double *) * shortperiods);	// *(pt_LGresults + i)
      for (int j = 0; j < shortperiods; j++)
	{
	  *(*(pt_LGresults + i) + j) = malloc (sizeof (double) * olsteps);	//*(*(pt_LGresults + i) + j)
	}
    }
  return pt_LGresults;
}

#ifdef DI
uint8_t **
allocatedidecisionarray (int nchannels, int nshortperiods)
{
  uint8_t **pt_shorttermdidecisionarray =
    malloc (sizeof (uint8_t *) * nchannels);
  for (int i = 0; i < nchannels; i++)
    {
      *(pt_shorttermdidecisionarray + i) =
	malloc (sizeof (uint8_t) * nshortperiods);
      //array could be longer than real audio so valgrind would blame this... initializing all to zero here is
      //possible but computationally intensive
    }
  //pt_shorttermdidecisionarray = NULL;
  return pt_shorttermdidecisionarray;
}

#endif

double **
allocatesc_shorttermarray (int nchannels, int nshortperiods)
{
  double **sc_shorttermarray = malloc (sizeof (double *) * nchannels);
  for (int i = 0; i < nchannels; i++)
    {
      *(sc_shorttermarray + i) = malloc (sizeof (double) * nshortperiods);
    }
  //pt_shorttermdidecisionarray = NULL;
  return sc_shorttermarray;
}


LG_Buf *
allocateLGBuffer (int samplenumber)
{
  LG_Buf *pt_LG_Buf;
  pt_LG_Buf = malloc (sizeof (LG_Buf));
  pt_LG_Buf->bufferA = malloc (sizeof (double) * samplenumber);
  pt_LG_Buf->bufferB = malloc (sizeof (double) * samplenumber);
  pt_LG_Buf->bufferLG = malloc (sizeof (double) * samplenumber);
  pt_LG_Buf->bufferSwap = malloc (sizeof (double) * samplenumber);
#ifdef DEBUG
  debuginit (pt_LG_Buf->bufferA, 4444.4444, samplenumber);
  debuginit (pt_LG_Buf->bufferB, 4444.4444, samplenumber);
  debuginit (pt_LG_Buf->bufferLG, 44444.4444, samplenumber);
  debuginit (pt_LG_Buf->bufferSwap, 4444.4444, samplenumber);
#endif
  return pt_LG_Buf;
}

void
debuginit (double *pt_toarray, double value, int arraylength)
{
  for (int i = 0; i < arraylength; i++)
    {
      pt_toarray[i] = value;
    }

}

int
freeLGBuffer (LG_Buf * pt_LG_Buf)
{
  free (pt_LG_Buf->bufferA);
  pt_LG_Buf->bufferA = NULL;
  free (pt_LG_Buf->bufferB);
  pt_LG_Buf->bufferB = NULL;
  free (pt_LG_Buf->bufferLG);
  pt_LG_Buf->bufferLG = NULL;
  free (pt_LG_Buf->bufferSwap);
  pt_LG_Buf->bufferSwap = NULL;
  free (pt_LG_Buf);
  pt_LG_Buf = NULL;
}

int
freeLGresultarray (int nchannels, int shortperiods, double ***lg_results)
{
  for (int i = 0; i < nchannels; i++)
    {
      for (int j = 0; j < shortperiods; j++)
	{
	  free (lg_results[i][j]);
	  lg_results[i][j] = NULL;
	}
      free (lg_results[i]);
      lg_results[i] = NULL;
    }
  free (lg_results);
  lg_results = NULL;
}


#ifdef DI

int
freedidecisionarray (int nchannels, uint8_t ** pt_shorttermdidecisionarray)
{
  for (int i = 0; i < nchannels; i++)
    {
      free (pt_shorttermdidecisionarray[i]);
      pt_shorttermdidecisionarray[i] = NULL;
    }
  free (pt_shorttermdidecisionarray);
  pt_shorttermdidecisionarray = NULL;
}
#endif

int
freesc_shorttermarray (int nchannels, double **sc_shorttermarray)
{
  for (int i = 0; i < nchannels; i++)
    {
      free (sc_shorttermarray[i]);
      sc_shorttermarray[i] = NULL;
    }
  free (sc_shorttermarray);
  sc_shorttermarray = NULL;
}




unsigned int
deintsamples (float *deintbuffer, int buffersizeframes, double *intbuffer,
	      int channel, int tot_channel_numb)
{
  unsigned int copied_samples = 0;
  int buffersizesamples = buffersizeframes * tot_channel_numb;

  /*
     This will copy buffersizesmpls of channel channel from intbuffer into deintbuffer.
   */
  for (int n = channel, m = 0; n < buffersizesamples;
       n += tot_channel_numb, m++)
    {
/*
#ifdef DEBUG

  printf("%.16f\n", intbuffer[n]);

#endif
*/
      deintbuffer[m] = (float) intbuffer[n];
      copied_samples = m; //this should be +1 and it should be possible to avoid it
/*      
#ifdef DEBUG
      printf(" %.8f\n", deintbuffer[m]);
#endif
  */    

    }


  return copied_samples + 1; //because index start with 0

}


#ifdef DI

						/* Conflate this with function with callback */

void *
di_worker_function (void *argstruct)
{

  struct WorkerArgs *thisWorkerArgs = (struct WorkerArgs *) argstruct;
  int output;
  if (thisWorkerArgs->is_eof_array[thisWorkerArgs->channel] == 0)
    {
      thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] =
	thisWorkerArgs->nsamples;
      if (thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] !=
	  thisWorkerArgs->fullbuffer_nsamples)
	{
	  thisWorkerArgs->is_eof_array[thisWorkerArgs->channel] = 1;
	}
    }
  else
    {
      thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] =
	(thisWorkerArgs->buffered_samples_array[thisWorkerArgs->channel] <
	 thisWorkerArgs->fullbuffer_nsamples) ?
	thisWorkerArgs->buffered_samples_array[thisWorkerArgs->channel] :
	thisWorkerArgs->fullbuffer_nsamples;
      thisWorkerArgs->buffered_samples_array[thisWorkerArgs->channel] -=
	thisWorkerArgs->samples_read_array[thisWorkerArgs->channel];
      memset (thisWorkerArgs->di_argbuffer, 0,
	      thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] *
	      sizeof (float));
    }
  if (thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] > 0)
    {
      output =
	di_process (thisWorkerArgs->di_array[thisWorkerArgs->channel],
		    thisWorkerArgs->di_argbuffer,
		    thisWorkerArgs->samples_read_array[thisWorkerArgs->
						       channel]);
      if (output == INVALID)
	{
	  thisWorkerArgs->buffered_samples_array[thisWorkerArgs->channel] +=
	    thisWorkerArgs->samples_read_array[thisWorkerArgs->channel];
	  /* consider non dialogue if output is invalid, otherwise values are not initialized */

	  thisWorkerArgs->shorttermarray_di[thisWorkerArgs->channel]
	    [thisWorkerArgs->shorttermindex] = (uint8_t) 0x00;


	}
      else
	{
	  thisWorkerArgs->shorttermarray_di[thisWorkerArgs->channel]
	    [thisWorkerArgs->shorttermindex] =
	    (uint8_t) ((output == SPEECH) ? 0x01 : 0x00);
	}

    }				// if (thisWorkerArgs->samples_read_array[thisWorkerArgs->channel] > 0)

  /* Here I have to free all possible things */

  /* Still I need to fill buffer with zeros if it is not a certain dimension. But most probably not here. */

  /* Free the memory allocated for DI */
  //    free(p_di);

  return 0;
}

#endif

#ifdef DI
void
print_di (int nchannels, int nshorttermperiods, uint8_t ** di_decisionarray)
{

  for (int i = 0; i < nchannels; i++)
    {
      printf ("Channel %d:\n", i);
      printf ("Period #:\t\tDecision:\n");
      for (int l = 0; l < nshorttermperiods; l++)
	{
	  printf ("%d\t\t%d\n", l, (int) di_decisionarray[i][l]);
	}
    }
}

#endif

#ifdef DI
						/*
						   chgatearray  - could be manually configured on the command line or automatically according to a Speech Content Threashold as per p. 17 of DolbyDialog Documentation. It is specified per channel. 0 is no gate, 1 is level gate, 2 is dialogue intelligence. At present level gate is not yet implemented

						   sc_staa -> short term avaraged array
						   stdda -> short term dialogue decision array
						   nch -> number of channels 
						   stp -> number of short term period number 
						   chgatearray -> channel gate configuration array

						 */

void
dolbydifinalcomputation (double **sc_staa, uint8_t ** stdda, int nch,
			 int stpn, int chgateconfarray[], struct Sum *ptSum)
{
  double spaccumulator = 0.0;
  int accumulatedperiods = 0;


  for (int i = 0; i < stpn; i++)
    {
      double intraperiodaccumulator = 0.0;
      int chgateor = 0;
      for (int j = 0; j < nch; j++)
	{
	  if ((chgateconfarray[j] == 2) && (stdda[j][i] == 1))
	    {			//At first I thought I would add all channels even if they do not have dialogue if only one of the channels hosting dialogue has actually dialogue.
	      intraperiodaccumulator += sc_staa[j][i];
	      if (stdda[j][i] == 1)
		{
		  printf ("%.10f\n", sc_staa[j][i]);
		}
	    }
	  if ((chgateconfarray[j] == 2) && (stdda[j][i] == 1))
	    {			//That is if channel configured for DI and DI decision "Speech/Other" is "Dialog" 
	      //spaccumulator += sc_staa[j][i];
	      chgateor = 1;

	    }			//if
	}			//for channel

      /*

         Here I should also insert the ~38 Leq(M) gate, or before?
         Take LKFS as model...
       */
      if (chgateor == 1)
	{
	  spaccumulator += intraperiodaccumulator;	// no if needed as intraperiod will be 0 otherwise
	  accumulatedperiods++;
	  //Debug information:
	  printf ("Accumulated a short period\t %.8f \t %.8f \n",
		  spaccumulator, intraperiodaccumulator);
	}
    }				//for short periods
  ptSum->dgleqm =
    20 * log10 (pow (spaccumulator / ((double) accumulatedperiods), 0.500)) +
    108.010299957;
  ptSum->dialoguepercentual =
    ((double) accumulatedperiods / (double) stpn) * 100;
}


#endif




void
lkfs_finalcomputation (LG * pt_lgctx, int *pt_chgateconf, int nchannels)
{
  int i_ch;
  //double gamma_A = leqmtosum(-70.0, 1);
  double gamma_A = -70.0;
  double gamma_r = -10.00;
  double gamma_R;

  /* expanding on chgateconf
     3 not included
     0 no gating
     1 level gating
     2 dialogue gating
   */

  double *ch_accumulator;	// first index channel, second threshold block


  ch_accumulator = malloc (sizeof (double) * pt_lgctx->totalsteps);	//or stepcounter ? or totalsteps? well totalsteps and shortperiod are calculated for a 5 hour feature in case of ffmpeg so use stepcounter later for the summing up




  int i_lgb = 0;		//level gate block index
  double LD_accum = 0.0;
  double LD_meas = 0.0;
  int gated_A_index = 0;
  int gated_R_index = 0;
  //int p_sp = 0; // short period index
  //int s_sd = 0; // subdivisions index
  //   for (int p_sp = 0; p_sp < pt_lgctx->shortperiods; p_sp++) {
  for (int p_sp = 0; p_sp < pt_lgctx->stepcounter / pt_lgctx->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx->subdivs; s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  /*
	     #ifdef DEBUG
	     printf ("Loudness (without Gating): %.4f\n", LD_meas);
	     #endif
	   */
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb++] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb++] = 0.0;
	    }			//
	}
    }

  //remainder
  if ((pt_lgctx->stepcounter % pt_lgctx->subdivs) != 0)
    {
      for (int s_sd = 0; s_sd < (pt_lgctx->stepcounter % pt_lgctx->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][pt_lgctx->stepcounter /
					      pt_lgctx->subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb++] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb++] = 0.0;
	    }			//
	}

    }

  // Relative Gating


  double gated_A_accum = 0.0;
  for (i_lgb = 0; i_lgb < pt_lgctx->stepcounter; i_lgb++)
    {
      gated_A_accum += ch_accumulator[i_lgb];
    }


  gamma_R =
    -0.691 + 10 * log10 (gated_A_accum / ((double) gated_A_index)) - 10;

  i_lgb = 0;

  for (int p_sp = 0; p_sp < pt_lgctx->stepcounter / pt_lgctx->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx->subdivs; s_sd++)
	{
	  //Relative gating
	  LD_accum = 0.0;
	  LD_meas = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb++] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb++] = 0.0;
	    }			//
	}
    }

  //remainder

  //remainder
  if ((pt_lgctx->stepcounter % pt_lgctx->subdivs) != 0)
    {
      for (int s_sd = 0; s_sd < (pt_lgctx->stepcounter % pt_lgctx->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0; 
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][pt_lgctx->stepcounter /
					      pt_lgctx->subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb++] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb++] = 0.0;
	    }			//
	}
    }

  double LKFS_accum = 0.0;
  double LKFS;

  for (i_lgb = 0; i_lgb < pt_lgctx->stepcounter; i_lgb++)
    {
      LKFS_accum += ch_accumulator[i_lgb];
    }
  LKFS = -0.691 + 10 * log10 (LKFS_accum / ((double) gated_R_index));

  printf ("LKFS: %.4f\n", LKFS);


  free (ch_accumulator);
}


#ifdef DI

void
lkfs_finalcomputation_withdolbydi (LG * pt_lgctx, int *pt_chgateconf,
				   int nchannels,
				   uint8_t ** stdda, double adlgthreshold)
{
  int i_ch;
  //double gamma_A = leqmtosum(-70.0, 1);
  double gamma_A = -70.0;
  double gamma_r = -10.00;
  double gamma_R;

  /* expanding on chgateconf
     3 not included
     0 no gating
     1 level gating
     2 dialogue gating
   */

  double *ch_accumulator;	// first index channel, second threshold block
  double *dich_accumulator;

  ch_accumulator = malloc (sizeof (double) * pt_lgctx->totalsteps);	//or stepcounter ? or totalsteps? well totalsteps and shortperiod are calculated for a 5 hour feature in case of ffmpeg so use stepcounter later for the summing up
  dich_accumulator = malloc (sizeof (double) * pt_lgctx->totalsteps);


  int digatedsignal = 0;
  int digatedcounter = 0;
  int digatedcounter_percent = 0;
  int i_lgb = 0;		//level gate block index
  double LD_accum = 0.0;
  double LD_meas = 0.0;
  double DI_accum = 0.0;
  double DI_meas = 0.0;
  int gated_A_index = 0;
  int gated_R_index = 0;
  //int p_sp = 0; // short period index
  //int s_sd = 0; // subdivisions index
  //   for (int p_sp = 0; p_sp < pt_lgctx->shortperiods; p_sp++) {
  for (int p_sp = 0; p_sp < pt_lgctx->stepcounter / pt_lgctx->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx->subdivs; s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  DI_accum = 0.0;

	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      if ((stdda[i_ch][p_sp] == 1) && (pt_chgateconf[i_ch] == 3))
		{
		  DI_accum += pt_lgctx->chgainconf[i_ch] * pt_lgctx->LGresultarray[i_ch][p_sp][s_sd];	// here chgain is not
		  digatedsignal = 1;
		}
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }

	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  DI_meas = -0.691 + 10 * log10 (DI_accum);
	  /*
	     #ifdef DEBUG
	     printf
	     ("Block Number: %08d, Loudness (without -70 Gating): %.4f\tDI Loudness (w/o -70 Gating): %.4f\n",
	     i_lgb, LD_meas, DI_meas);
	     #endif
	   */
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }
	  if (DI_meas > gamma_A)
	    {
	      dich_accumulator[i_lgb] = DI_accum;
	      digatedcounter++; // In this implementation counter is increased after absolute level gating, this differs from p. 17. No re-reading
	      // I think that according to the integration documentation one has to have two counters, one for the mean the other for the percentage
	      // because fpr the mean the counter must be incremented here, otherwise the measure get down arbitrarily for the blocks below absolute threshold
	    }
	  else
	    {
	      dich_accumulator[i_lgb] = 0.0;
	    }
	  i_lgb++;
	  if (digatedsignal == 1)
	    {
	      digatedcounter_percent++;	// This is the same as in p.17 Dolby Dialogue Intelligence Reference Code User's Guide
	      digatedsignal = 0;
	    }
	}
    }

  //remainder
  if ((pt_lgctx->stepcounter % pt_lgctx->subdivs) != 0)
    {
      for (int s_sd = 0; s_sd < (pt_lgctx->stepcounter % pt_lgctx->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0;
	  DI_accum = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      if ((stdda[i_ch][pt_lgctx->stepcounter / pt_lgctx->subdivs] ==
		   1) && (pt_chgateconf[i_ch] == 3))
		{
		  DI_accum += pt_lgctx->chgainconf[i_ch] * pt_lgctx->LGresultarray[i_ch][pt_lgctx->stepcounter / pt_lgctx->subdivs][s_sd];	// here chgain is not
		  digatedsignal = 1;
		}

	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][pt_lgctx->stepcounter /
					      pt_lgctx->subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  DI_meas = -0.691 + 10 * log10 (DI_accum);
	  /*
	     #ifdef DEBUG
	     printf
	     ("Block Number: %08d, Loudness (without -70 Gating): %.4f\tDI Loudness (w/o -70 Gating): %.4f\n",
	     i_lgb, LD_meas, DI_meas);
	     #endif
	   */
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  if (DI_meas > gamma_A)
	    {
	      dich_accumulator[i_lgb] = DI_accum;
	      digatedcounter++; // In this implementation counter is increased after absolute level gating, this differs from p. 17. But I changed my mind rereading and added a second counter see above.
	    }
	  else
	    {
	      dich_accumulator[i_lgb] = 0.0;
	    }
	  i_lgb++;
	  if (digatedsignal == 1)
	    {
	      digatedcounter_percent++;	// This is the same as in p.17 Dolby Dialogue Intelligence Reference Code User's Guide
	      digatedsignal = 0;
	    }

	}

    }

  // Relative Gating


  double gated_A_accum = 0.0;
  for (i_lgb = 0; i_lgb < pt_lgctx->stepcounter; i_lgb++)
    {
      gated_A_accum += ch_accumulator[i_lgb];
    }


  gamma_R =
    -0.691 + 10 * log10 (gated_A_accum / ((double) gated_A_index)) - 10;

  i_lgb = 0;

  for (int p_sp = 0; p_sp < pt_lgctx->stepcounter / pt_lgctx->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx->subdivs; s_sd++)
	{
	  //Relative gating
	  LD_accum = 0.0;
	  LD_meas = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  i_lgb++;
	}
    }

  //remainder

  //remainder
  if ((pt_lgctx->stepcounter % pt_lgctx->subdivs) != 0)
    {
      for (int s_sd = 0; s_sd < (pt_lgctx->stepcounter % pt_lgctx->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0; 
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx->chgainconf[i_ch] *
		pt_lgctx->LGresultarray[i_ch][pt_lgctx->stepcounter /
					      pt_lgctx->subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = -0.691 + 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  i_lgb++;
	}
    }

  double DI_LKFS_accum = 0.0;
  double LKFS_accum = 0.0;
  double LKFS;
  double DILKFS;
  double dialoguepercentage;

  for (i_lgb = 0; i_lgb < pt_lgctx->stepcounter; i_lgb++)
    {
      LKFS_accum += ch_accumulator[i_lgb];
    }
  LKFS = -0.691 + 10 * log10 (LKFS_accum / ((double) gated_R_index));

  printf ("LKFS: %.4f\n", LKFS);

  dialoguepercentage =
    ((double) digatedcounter_percent) / ((double) pt_lgctx->stepcounter) * 100.00;

  printf ("Speech Percentage is: %.2f %%\n", dialoguepercentage);

  if (dialoguepercentage >= adlgthreshold)
    {
      for (i_lgb = 0; i_lgb < pt_lgctx->stepcounter; i_lgb++)
	{
	  DI_LKFS_accum += dich_accumulator[i_lgb];
	}
      DILKFS =
	-0.691 + 10 * log10 (DI_LKFS_accum / ((double) digatedcounter));
      printf ("LKFS(DI): %.2f\n", DILKFS);
    }
  free (ch_accumulator);
  free (dich_accumulator);
}


#endif


#ifdef DI

void
dolbydifinalcomputation2 (LGLeqM * pt_lgctx_leqmdi, int *pt_chgateconf,
			  int nchannels,
			  uint8_t ** stdda, double adlgthreshold)
{

  int i_ch;
  //double gamma_A = leqmtosum(-70.0, 1);
  double gamma_A = -70.0;
  double gamma_r = -10.00;
  double gamma_R;

  /* expanding on chgateconf
     3 not included
     0 no gating
     1 level gating
     2 dialogue gating
   */

  double *ch_accumulator;	// first index channel, second threshold block
  double *dich_accumulator;

  ch_accumulator = malloc (sizeof (double) * pt_lgctx_leqmdi->totalsteps);	//or stepcounter ? or totalsteps? well totalsteps and shortperiod are calculated for a 5 hour feature in case of ffmpeg so use stepcounter later for the summing up
  dich_accumulator = malloc (sizeof (double) * pt_lgctx_leqmdi->totalsteps);


  int digatedsignal = 0;
  int digatedcounter = 0;
  int digatedcounter_percent = 0;
  int i_lgb = 0;		//level gate block index
  double LD_accum = 0.0;
  double LD_meas = 0.0;
  double DI_accum = 0.0;
  double DI_meas = 0.0;
  int gated_A_index = 0;
  int gated_R_index = 0;
  //int p_sp = 0; // short period index
  //int s_sd = 0; // subdivisions index
  //   for (int p_sp = 0; p_sp < pt_lgctx->shortperiods; p_sp++) {
  for (int p_sp = 0;
       p_sp < pt_lgctx_leqmdi->stepcounter / pt_lgctx_leqmdi->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx_leqmdi->subdivs; s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  DI_accum = 0.0;

	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      if ((stdda[i_ch][p_sp] == 1) && (pt_chgateconf[i_ch] == 3))
		{
		  DI_accum += pt_lgctx_leqmdi->chgainconf[i_ch] * pt_lgctx_leqmdi->LGresultarray[i_ch][p_sp][s_sd];	// here chgain is not
		  digatedsignal = 1;
		}
	      LD_accum +=
		pt_lgctx_leqmdi->chgainconf[i_ch] *
		pt_lgctx_leqmdi->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }

	  //Loudness without Gating (single block)
	  LD_meas = 10 * log10 (LD_accum);
	  DI_meas = 10 * log10 (DI_accum);
	  /*
	     #ifdef DEBUG
	     printf
	     ("Block Number: %08d, Loudness (without -70 Gating): %.4f\tDI Loudness (w/o -70 Gating): %.4f\n",
	     i_lgb, LD_meas, DI_meas);
	     #endif
	   */
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }
	  if (DI_meas > gamma_A)
	    {
	      dich_accumulator[i_lgb] = DI_accum;
	      digatedcounter++; // In this implementation counter is increased after absolute level gating, this differs from p. 17
	    }
	  else
	    {
	      dich_accumulator[i_lgb] = 0.0;
	    }
	  i_lgb++;
	  if (digatedsignal == 1)
	    {
	      digatedcounter_percent++;	// This is the same as in p.17 Dolby Dialogue Intelligence Reference Code User's Guide
	      digatedsignal = 0;
	    }
	}
    }

  //remainder
  if ((pt_lgctx_leqmdi->stepcounter % pt_lgctx_leqmdi->subdivs) != 0)
    {
      for (int s_sd = 0;
	   s_sd < (pt_lgctx_leqmdi->stepcounter % pt_lgctx_leqmdi->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0;
	  DI_accum = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      if ((stdda[i_ch]
		   [pt_lgctx_leqmdi->stepcounter /
		    pt_lgctx_leqmdi->subdivs] == 1)
		  && (pt_chgateconf[i_ch] == 3))
		{
		  DI_accum += pt_lgctx_leqmdi->chgainconf[i_ch] * pt_lgctx_leqmdi->LGresultarray[i_ch][pt_lgctx_leqmdi->stepcounter / pt_lgctx_leqmdi->subdivs][s_sd];	// here chgain is not
		  digatedsignal = 1;
		}

	      LD_accum +=
		pt_lgctx_leqmdi->chgainconf[i_ch] *
		pt_lgctx_leqmdi->LGresultarray[i_ch][pt_lgctx_leqmdi->
						     stepcounter /
						     pt_lgctx_leqmdi->
						     subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = 10 * log10 (LD_accum);
	  DI_meas = 10 * log10 (DI_accum);
	  /*
	     #ifdef DEBUG
	     printf
	     ("Block Number: %08d, Loudness (without -70 Gating): %.4f\tDI Loudness (w/o -70 Gating): %.4f\n",
	     i_lgb, LD_meas, DI_meas);
	     #endif
	   */
	  if (LD_meas > gamma_A)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_A_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  if (DI_meas > gamma_A)
	    {
	      dich_accumulator[i_lgb] = DI_accum;
	      digatedcounter++; // In this implementation counter is increased after absolute level gating, this differs from p. 17
	    }
	  else
	    {
	      dich_accumulator[i_lgb] = 0.0;
	    }
	  i_lgb++;
	  if (digatedsignal == 1)
	    {
	      digatedcounter_percent++;	// This is the same as in p.17 Dolby Dialogue Intelligence Reference Code User's Guide
	      digatedsignal = 0;
	    }

	}

    }

  // Relative Gating


  double gated_A_accum = 0.0;
  for (i_lgb = 0; i_lgb < pt_lgctx_leqmdi->stepcounter; i_lgb++)
    {
      gated_A_accum += ch_accumulator[i_lgb];
    }


  gamma_R =
    10 * log10 (gated_A_accum / ((double) gated_A_index)) - 10;

  i_lgb = 0;

  for (int p_sp = 0;
       p_sp < pt_lgctx_leqmdi->stepcounter / pt_lgctx_leqmdi->subdivs; p_sp++)
    {
      for (int s_sd = 0; s_sd < pt_lgctx_leqmdi->subdivs; s_sd++)
	{
	  //Relative gating
	  LD_accum = 0.0;
	  LD_meas = 0.0;
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx_leqmdi->chgainconf[i_ch] *
		pt_lgctx_leqmdi->LGresultarray[i_ch][p_sp][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  i_lgb++;
	}
    }

  //remainder

  //remainder
  if ((pt_lgctx_leqmdi->stepcounter % pt_lgctx_leqmdi->subdivs) != 0)
    {
      for (int s_sd = 0;
	   s_sd < (pt_lgctx_leqmdi->stepcounter % pt_lgctx_leqmdi->subdivs);
	   s_sd++)
	{
	  //Absolute gating gating
	  LD_accum = 0.0;
	  //LD_meas = 0.0; 
	  for (i_ch = 0; i_ch < nchannels; i_ch++)
	    {
	      LD_accum +=
		pt_lgctx_leqmdi->chgainconf[i_ch] *
		pt_lgctx_leqmdi->LGresultarray[i_ch][pt_lgctx_leqmdi->
						     stepcounter /
						     pt_lgctx_leqmdi->
						     subdivs][s_sd];

	      // see p. 5 of ITU 1770-4
	    }
	  //Loudness without Gating (single block)
	  LD_meas = 10 * log10 (LD_accum);
	  if (LD_meas > gamma_R)
	    {
	      ch_accumulator[i_lgb] = LD_accum;
	      gated_R_index++;
	    }
	  else
	    {
	      ch_accumulator[i_lgb] = 0.0;
	    }			//
	  i_lgb++;
	}
    }

  double DI_LGLEQM_accum = 0.0;
  double LGLEQM_accum = 0.0;
  double LGLEQM;
  double DI_LGLEQM;
  double dialoguepercentage;

  for (i_lgb = 0; i_lgb < pt_lgctx_leqmdi->stepcounter; i_lgb++)
    {
      LGLEQM_accum += ch_accumulator[i_lgb];
    }
  LGLEQM = 10 * log10 (LGLEQM_accum / ((double) gated_R_index));	// not taking the square root because multiplying by 10, it is indeed power

  printf ("Leq(M,LG)FS: %.4f\n", LGLEQM);
  printf ("Leq(M,LG): %.4f\n", LGLEQM + 108.010299957);
  dialoguepercentage =
    ((double) digatedcounter_percent) / ((double) pt_lgctx_leqmdi->stepcounter) *
    100.00;

  printf ("Speech Percentage: %.2f %%\n", dialoguepercentage);

  if (dialoguepercentage >= adlgthreshold)
    {
      for (i_lgb = 0; i_lgb < pt_lgctx_leqmdi->stepcounter; i_lgb++)
	{
	  DI_LGLEQM_accum += dich_accumulator[i_lgb];
	}
      DI_LGLEQM = 10 * log10 (DI_LGLEQM_accum / ((double) digatedcounter));	// it is power
      printf ("Leq(M,DI)FS: %.4f\n", DI_LGLEQM);
      printf ("Leq(M,DI): %.4f\n", DI_LGLEQM + 108.010299957);
    }
  free (ch_accumulator);
  free (dich_accumulator);



}


#endif


void
levelgatefinalcomputation (double **sc_staa, double linearthreshold, int nch,
			   int stpn, struct Sum *ptSum)
{

  /* Threshold is not in dB, but already transformed for comparison */

  double spaccumulator = 0.0;
  int accumulatedperiods = 0;


  for (int i = 0; i < stpn; i++)
    {
      double intraperiodaccumulator = 0.0;
      int chgateor = 0;
      for (int j = 0; j < nch; j++)
	{
	  intraperiodaccumulator += sc_staa[j][i];
	}
      if (intraperiodaccumulator >= linearthreshold)
	{
	  printf ("Added to simple level gate accumulator\n");
	  spaccumulator += intraperiodaccumulator;
	  accumulatedperiods++;
	}
    }


  ptSum->lgleqm =
    20 * log10 (pow (spaccumulator / ((double) accumulatedperiods), 0.500)) +
    108.010299957;


}


void
adaptivegateselection (void)
{

}


double
leqmtosum (double leqm, double ref)
{
  /* this is used to transfom the argument of --levelgate into a linear level ready for comparison in the accumulator 
     ref would be most of the time 20uPa -- and here again the trouble with 20uPa...
   */

  return powf (10, (leqm / 20)) * ref;

}


								/* The following comes from Rec. ITU-R BS. 1770-4 pp.4-5 */


int
K_filter_stage1 (double *smp_out, double *smp_in, int nsamples,
		 coeff * coeff_ctx)
{
  /* I have to limit it for the first 2 samples */
  /* This is coming from ITU ... */
  /* This coefficient are  for 48kHz sample rate */

  for (int i = 0; i < nsamples; i++)
    {
      switch (i)
	{
	case 0:
	  smp_out[i] = coeff_ctx->hsb0 * smp_in[i];
	  break;
	case 1:
	  smp_out[i] =
	    coeff_ctx->hsb0 * smp_in[i] + coeff_ctx->hsb1 * smp_in[i - 1] -
	    coeff_ctx->hsa1 * smp_out[i - 1];
	  break;
	default:
	  smp_out[i] =
	    coeff_ctx->hsb0 * smp_in[i] + coeff_ctx->hsb1 * smp_in[i - 1] +
	    coeff_ctx->hsb2 * smp_in[i - 2] - coeff_ctx->hsa1 * smp_out[i -
									1] -
	    coeff_ctx->hsa2 * smp_out[i - 2];
	  break;
	}
    }
}



int
K_filter_stage2 (double *smp_out, double *smp_in, int nsamples,
		 coeff * coeff_ctx)
{
  /* I have to limit it for the first 2 samples */
  /* This coefficient are only for 48kHz sample rate */
  for (int i = 0; i < nsamples; i++)
    {
      switch (i)
	{
	case 0:
	  smp_out[i] = coeff_ctx->hpb0 * smp_in[i];
	  break;
	case 1:
	  smp_out[i] =
	    coeff_ctx->hpb0 * smp_in[i] + coeff_ctx->hpb1 * smp_in[i - 1] -
	    coeff_ctx->hpa1 * smp_out[i - 1];
	  break;
	default:
	  smp_out[i] =
	    coeff_ctx->hpb0 * smp_in[i] + coeff_ctx->hpb1 * smp_in[i - 1] +
	    coeff_ctx->hpb2 * smp_in[i - 2] - coeff_ctx->hpa1 * smp_out[i -
									1] -
	    coeff_ctx->hpa2 * smp_out[i - 2];
	  break;
	}
    }
}







								/* Work on Gating Block for Leq(M) or LUFS */


int
calcSampleStepLG (float percentOverlap, int samplerate, int LGbufferms)
{
  int sampleStep =
    (int) ((((float) samplerate) * ((float) LGbufferms) / 1000.00) *
	   ((100 - percentOverlap) / 100.00));
#ifdef DEBUG
  printf ("Level Gate step in samples: %d.\n", sampleStep);
#endif
  return sampleStep;
}



								/*

								   %%%%%%
								   Order 5 FS 44100

								   A =

								   Columns 1 through 5:

								   1.0000000000000000  -1.5224995723629664   1.3617953870010380  -0.7794603877415162   0.2773974331876455

								   Column 6:

								   -0.0477648119172564

								   >> B
								   B =

								   Columns 1 through 5:

								   0.4034108659797224   0.0675046624145518  -0.3122917473135974  -0.1471391464872613  -0.0173711282192394

								   Column 6:

								   0.0101026340442429

								   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


								   New Oder 4  FS 48000

								   >> A
								   A =

								   1.000000000000000  -1.509001373787255   1.263289551410631  -0.603342410310913   0.143867967944599

								   >> B
								   B =

								   0.3183679456956268   0.1494237633315084  -0.2095809948670236  -0.1892237517582008  -0.0603259946368122

								   >>

								   New Oder 5 FS 48000

								   A =

								   Columns 1 through 5:

								   1.0000000000000000  -1.6391291074367320   1.5160386192837869  -0.8555167646249104   0.2870466545317107

								   Column 6:

								   -0.0428951718612053

								   >> B
								   B =

								   Columns 1 through 4:

								   0.31837346242469328   0.10800452155339044  -0.21106344349319428  -0.15438275853192485

								   Columns 5 and 6:

								   -0.05130596901975942  -0.00518224535906041

								   >>


								   %%%%%%%%%%%
								   oder 7 Fs = 96000

								   A =

								   Columns 1 through 5:

								   1.000000000000000  -4.065513964912803   7.729630359954706  -9.008501669250943   7.042491931110629

								   Columns 6 through 8:

								   -3.718680112548192   1.223024516229922  -0.190932048752515

								   >> B
								   B =

								   Columns 1 through 4:

								   2.20508445245658e-02   4.75013833045082e-03  -1.14732527362953e-02   9.84044708222169e-04

								   Columns 5 through 8:

								   -5.97112489049753e-03  -8.65702130320201e-03  -2.10596182798782e-03   4.28003720960765e-04

								   >>

								   Order 9 FS 192000

								   A =

								   Columns 1 through 5:

								   1.000000000000000   -6.637252076760877   19.972276555714540  -35.816564526532424   42.281935151372885

								   Columns 6 through 10:

								   -34.193652600767841   19.040936622298062   -7.090648293919744    1.616538001967351   -0.173496773056347

								   >> B
								   B =

								   Columns 1 through 4:

								   1.44318321996676e-04   2.22020582676945e-04  -1.04877188873611e-04  -3.01386180177350e-04

								   Columns 5 through 8:

								   -7.03016165414926e-06   2.04514103920991e-04   5.85534078264994e-05  -9.99320871937674e-05

								   Columns 9 and 10:

								   -9.02747116074520e-05  -2.72118221944114e-05


								 */




int
M_filter (double *smp_out, double *smp_in, int samples, int samplerate)
{
  if (samplerate == 44100)
    {
      for (int i = 0; i < samples; i++)
	{

	  // Order 5
	  switch (i)
	    {
	    case 0:
	      smp_out[i] = 0.4034108659797224 * smp_in[i];
	      break;
	    case 1:
	      smp_out[i] =
		0.4034108659797224 * smp_in[i] +
		0.0675046624145518 * smp_in[i - 1] +
		1.5224995723629664 * smp_out[i - 1];
	      break;
	    case 2:
	      smp_out[i] =
		0.4034108659797224 * smp_in[i] +
		0.0675046624145518 * smp_in[i - 1] -
		0.3122917473135974 * smp_in[i - 2] +
		1.5224995723629664 * smp_out[i - 1] -
		1.3617953870010380 * smp_out[i - 2];
	      break;
	    case 3:
	      smp_out[i] =
		0.4034108659797224 * smp_in[i] +
		0.0675046624145518 * smp_in[i - 1] -
		0.3122917473135974 * smp_in[i - 2] -
		0.1471391464872613 * smp_in[i - 3] +
		1.5224995723629664 * smp_out[i - 1] -
		1.3617953870010380 * smp_out[i - 2] +
		0.7794603877415162 * smp_out[i - 3];
	      break;
	    case 4:
	      smp_out[i] =
		0.4034108659797224 * smp_in[i] +
		0.0675046624145518 * smp_in[i - 1] -
		0.3122917473135974 * smp_in[i - 2] -
		0.1471391464872613 * smp_in[i - 3] -
		0.0173711282192394 * smp_in[i - 4] +
		1.5224995723629664 * smp_out[i - 1] -
		1.3617953870010380 * smp_out[i - 2] +
		0.7794603877415162 * smp_out[i - 3] -
		0.2773974331876455 * smp_out[i - 4];
	      break;
	    default:
	      smp_out[i] =
		0.4034108659797224 * smp_in[i] +
		0.0675046624145518 * smp_in[i - 1] -
		0.3122917473135974 * smp_in[i - 2] -
		0.1471391464872613 * smp_in[i - 3] -
		0.0173711282192394 * smp_in[i - 4] +
		0.0101026340442429 * smp_in[i - 5] +
		1.5224995723629664 * smp_out[i - 1] -
		1.3617953870010380 * smp_out[i - 2] +
		0.7794603877415162 * smp_out[i - 3] -
		0.2773974331876455 * smp_out[i - 4] +
		0.0477648119172564 * smp_out[i - 5];
	      break;
	    }



	}

    }
  else if (samplerate == 48000)
    {				// if (samplerate == 44100)

      /* This coefficient are only for 48kHz sample rate */
      for (int i = 0; i < samples; i++)
	{
	  /* // Order 4
	     switch(i) {
	     case 0:
	     smp_out[i] =  0.3160057437484637  * smp_in[i];
	     break;
	     case 1:
	     smp_out[i] = 0.3183679456956268 * smp_in[i] +  0.1494237633315084 * smp_in[i-1] + 1.509001373787255 * smp_out[i-1];
	     break;
	     case 2:
	     smp_out[i] = 0.3183679456956268 * smp_in[i] +  0.1494237633315084 * smp_in[i-1] - 0.2095809948670236 * smp_in[i-2] + 1.509001373787255 * smp_out[i-1] - 1.263289551410631 * smp_out[i-2];
	     break;
	     case 3:
	     smp_out[i] = 0.3183679456956268 * smp_in[i] +  0.1494237633315084 * smp_in[i-1] - 0.2095809948670236 * smp_in[i-2] - 0.1892237517582008 * smp_in[i-3] + 1.509001373787255 * smp_out[i-1] - 1.263289551410631 * smp_out[i-2] + 0.603342410310913  * smp_out[i-3];
	     break;
	     default:
	     smp_out[i] = 0.3183679456956268 * smp_in[i] +  0.1494237633315084 * smp_in[i-1] - 0.2095809948670236 * smp_in[i-2] - 0.1892237517582008 * smp_in[i-3]  -0.0603259946368122 * smp_in[i-4] + 1.509001373787255 * smp_out[i-1] - 1.263289551410631 * smp_out[i-2] + 0.603342410310913  * smp_out[i-3] - 0.143867967944599 * smp_out[i-4];
	     break;
	     } */
	  // Order 5
	  switch (i)
	    {
	    case 0:
	      smp_out[i] = 0.31837346242469328 * smp_in[i];
	      break;
	    case 1:
	      smp_out[i] =
		0.31837346242469328 * smp_in[i] +
		0.10800452155339044 * smp_in[i - 1] +
		1.6391291074367320 * smp_out[i - 1];
	      break;
	    case 2:
	      smp_out[i] =
		0.31837346242469328 * smp_in[i] +
		0.10800452155339044 * smp_in[i - 1] -
		0.21106344349319428 * smp_in[i - 2] +
		1.6391291074367320 * smp_out[i - 1] -
		1.5160386192837869 * smp_out[i - 2];
	      break;
	    case 3:
	      smp_out[i] =
		0.31837346242469328 * smp_in[i] +
		0.10800452155339044 * smp_in[i - 1] -
		0.21106344349319428 * smp_in[i - 2] -
		0.15438275853192485 * smp_in[i - 3] +
		1.6391291074367320 * smp_out[i - 1] -
		1.5160386192837869 * smp_out[i - 2] +
		0.8555167646249104 * smp_out[i - 3];
	      break;
	    case 4:
	      smp_out[i] =
		0.31837346242469328 * smp_in[i] +
		0.10800452155339044 * smp_in[i - 1] -
		0.21106344349319428 * smp_in[i - 2] -
		0.15438275853192485 * smp_in[i - 3] -
		0.05130596901975942 * smp_in[i - 4] +
		1.6391291074367320 * smp_out[i - 1] -
		1.5160386192837869 * smp_out[i - 2] +
		0.8555167646249104 * smp_out[i - 3] -
		0.2870466545317107 * smp_out[i - 4];
	      break;
	    default:
	      smp_out[i] =
		0.31837346242469328 * smp_in[i] +
		0.10800452155339044 * smp_in[i - 1] -
		0.21106344349319428 * smp_in[i - 2] -
		0.15438275853192485 * smp_in[i - 3] -
		0.05130596901975942 * smp_in[i - 4] -
		0.00518224535906041 * smp_in[i - 5] +
		1.6391291074367320 * smp_out[i - 1] -
		1.5160386192837869 * smp_out[i - 2] +
		0.8555167646249104 * smp_out[i - 3] -
		0.2870466545317107 * smp_out[i - 4] +
		0.0428951718612053 * smp_out[i - 5];
	      break;
	    }



	}			// for
    }
  else if (samplerate == 96000)
    {				// if (samplerate == 48000)
      for (int i = 0; i < samples; i++)
	{
	  switch (i)
	    {
	    case 0:
	      smp_out[i] = 2.20508445245658E-2 * smp_in[i];
	      break;
	    case 1:
	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] +
		4.065513964912803 * smp_out[i - 1];
	      break;
	    case 2:
	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2];
	      break;
	    case 3:

	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		9.84044708222169E-4 * smp_in[i - 3] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2] +
		9.008501669250943 * smp_out[i - 3];
	      break;
	    case 4:
	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		9.84044708222169E-4 * smp_in[i - 3] -
		5.97112489049753E-3 * smp_in[i - 4] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2] +
		9.008501669250943 * smp_out[i - 3] -
		7.042491931110629 * smp_out[i - 4];
	      break;

	    case 5:

	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		9.84044708222169E-4 * smp_in[i - 3] -
		5.97112489049753E-3 * smp_in[i - 4] -
		8.65702130320201E-3 * smp_in[i - 5] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2] +
		9.008501669250943 * smp_out[i - 3] -
		7.042491931110629 * smp_out[i - 4] +
		3.718680112548192 * smp_out[i - 5];
	      break;

	    case 6:

	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		9.84044708222169E-4 * smp_in[i - 3] -
		5.97112489049753E-3 * smp_in[i - 4] -
		8.65702130320201E-3 * smp_in[i - 5] -
		2.10596182798782E-3 * smp_in[i - 6] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2] +
		9.008501669250943 * smp_out[i - 3] -
		7.042491931110629 * smp_out[i - 4] +
		3.718680112548192 * smp_out[i - 5] -
		1.223024516229922 * smp_out[i - 6];
	      break;

	    default:
	      smp_out[i] =
		2.20508445245658E-2 * smp_in[i] +
		4.75013833045082E-3 * smp_in[i - 1] -
		1.14732527362953E-2 * smp_in[i - 2] +
		9.84044708222169E-4 * smp_in[i - 3] -
		5.97112489049753E-3 * smp_in[i - 4] -
		8.65702130320201E-3 * smp_in[i - 5] -
		2.10596182798782E-3 * smp_in[i - 6] +
		4.28003720960765E-4 * smp_in[i - 7] +
		4.065513964912803 * smp_out[i - 1] -
		7.729630359954706 * smp_out[i - 2] +
		9.008501669250943 * smp_out[i - 3] -
		7.042491931110629 * smp_out[i - 4] +
		3.718680112548192 * smp_out[i - 5] -
		1.223024516229922 * smp_out[i - 6] +
		0.190932048752515 * smp_out[i - 7];
	      break;

	    }
	}
    }
  else if (samplerate == 192000)
    {				// if (samplerate == 96000)
      for (int i = 0; i < samples; i++)
	{
	  switch (i)
	    {
	    case 0:
	      smp_out[i] = 1.44318321996676E-4 * smp_in[i];
	      break;
	    case 1:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] +
		6.637252076760877 * smp_out[i - 1];
	      break;
	    case 2:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2];
	      break;
	    case 3:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3];
	      break;
	    case 4:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4];
	      break;

	    case 5:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		2.04514103920991E-4 * smp_in[i - 5] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4] +
		34.193652600767841 * smp_out[i - 5];
	      break;

	    case 6:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		2.04514103920991E-4 * smp_in[i - 5] +
		5.85534078264994E-5 * smp_in[i - 6] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4] +
		34.193652600767841 * smp_out[i - 5] -
		19.040936622298062 * smp_out[i - 6];
	      break;
	    case 7:

	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		2.04514103920991E-4 * smp_in[i - 5] +
		5.85534078264994E-5 * smp_in[i - 6] -
		9.99320871937674E-5 * smp_in[i - 7] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4] +
		34.193652600767841 * smp_out[i - 5] -
		19.040936622298062 * smp_out[i - 6] +
		7.090648293919744 * smp_out[i - 7];
	      break;
	    case 8:

	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		2.04514103920991E-4 * smp_in[i - 5] +
		5.85534078264994E-5 * smp_in[i - 6] -
		9.99320871937674E-5 * smp_in[i - 7] -
		9.02747116074520E-5 * smp_in[i - 8] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4] +
		34.193652600767841 * smp_out[i - 5] -
		19.040936622298062 * smp_out[i - 6] +
		7.090648293919744 * smp_out[i - 7] -
		1.616538001967351 * smp_out[i - 8];
	      break;

	    default:
	      smp_out[i] =
		1.44318321996676E-4 * smp_in[i] +
		2.22020582676945E-4 * smp_in[i - 1] -
		1.04877188873611E-4 * smp_in[i - 2] -
		3.01386180177350E-4 * smp_in[i - 3] -
		7.03016165414926E-6 * smp_in[i - 4] +
		2.04514103920991E-4 * smp_in[i - 5] +
		5.85534078264994E-5 * smp_in[i - 6] -
		9.99320871937674E-5 * smp_in[i - 7] -
		9.02747116074520E-5 * smp_in[i - 8] -
		2.72118221944114E-5 * smp_in[i - 9] +
		6.637252076760877 * smp_out[i - 1] -
		19.972276555714540 * smp_out[i - 2] +
		35.816564526532424 * smp_out[i - 3] -
		42.281935151372885 * smp_out[i - 4] +
		34.193652600767841 * smp_out[i - 5] -
		19.040936622298062 * smp_out[i - 6] +
		7.090648293919744 * smp_out[i - 7] -
		1.616538001967351 * smp_out[i - 8] +
		0.173496773056347 * smp_out[i - 9];
	      break;

	    }

	}
    }
}


TruePeak *
init_truepeak_ctx (int ch, int os, int taps)
{
  TruePeak *tp = malloc (sizeof (TruePeak));
  tp->vector = malloc (sizeof (double) * ch);
  for (int i = 0; i < ch; i++)
    {
      tp->vector[i] = 0.0;
    }
  tp->oversampling_ratio = os;
  tp->filtertaps = taps;
  //Coefficients will be initialized by a dedicated function
  return tp;
}

int
freetruepeak (TruePeak * tp)
{
  free (tp->vector);
  tp->vector = NULL;
  free (tp->filter_coeffs);
  tp->filter_coeffs = NULL;
  free (tp);
}


double
truepeakcheck (double *in_buf, int ns, double truepeak, int os_ratio,
	       int filtertaps, double *coeff_vector)
{

  /*


     Samples are from a single channel

     - scale by -12.04 dB
     - oversample x4
     - lowpass filter
     - rectify
     - 20 log
     - scale by +12.04

   */


  double *os_buffer = malloc (sizeof (double) * ns * os_ratio);
  double *interp_buffer = calloc (ns * os_ratio, sizeof (double));
  double att_in_lin = pow (10, (-12.04 / 20));
  //double boost_out = +12.04;
  /*
     double poly_phase1 [] = { 0.0017089843750, 0.0109863281250, −0.0196533203125, 0.0332031250000, −0.0594482421875, 0.1373291015625,0.9721679687500,  −0.1022949218750, 0.0476074218750, −0.0266113281250, 0.0148925781250, −0.0083007812500};
     double poly_phase2 [] = {−0.0291748046875,0.0292968750000, −0.0517578125000, 0.0891113281250, −0.1665039062500,  0.4650878906250, 0.7797851562500,   −0.2003173828125, 0.1015625000000, −0.0582275390625, 0.0330810546875, −0.0189208984375};
     double poly_phase3 [] = {−0.0189208984375,   0.0330810546875, −0.0582275390625,  0.1015625000000,  −0.2003173828125,  0.7797851562500,  0.4650878906250,  −0.1665039062500, 0.0891113281250,  −0.0517578125000, 0.0292968750000,  −0.0291748046875};
     double poly_phase4 [] = { −0.0083007812500, 0.0148925781250, −0.0266113281250, 0.0476074218750, −0.1022949218750, 0.9721679687500, 0.1373291015625, −0.0594482421875, 0.0332031250000, −0.0196533203125, 0.0109863281250, 0.00170898437};

   */

  /*attenuate and oversample (extend filling with zeros) buffer */
  int i;
  int j;
  for (i = 0, j = 0; i < ns * os_ratio; i += os_ratio, j++)
    {
      os_buffer[i] = att_in_lin * in_buf[j];
      for (int k = 1; k < os_ratio; k++)
	{
	  os_buffer[i + k] = 0.0;
	}

    }

  /* low pass filter */

  for (int i = 0; i < (ns * os_ratio); i++)
    {
      for (int j = 0; j < min (i + 1, filtertaps); j++)
	{
	  interp_buffer[i] += coeff_vector[j] * os_buffer[i - j];
	  /*
	     #ifdef DEBUG
	     printf("Sample: %08d\tcoeff.: %03d\t value: %.10f\n", i, j, os_buffer[i]);
	     #endif
	   */
	}
    }

  for (int m = 0; m < (ns * os_ratio); m++)
    {
      /*
         #ifdef DEBUG
         printf("Old Truepeak: %.4f\n", truepeak);
         printf("Sample value: %.8f\n", os_buffer[m]);
         #endif       
       */
      truepeak = max (interp_buffer[m] * interp_buffer[m], truepeak);	//rectify

    }

  free (os_buffer);
  free (interp_buffer);
  return truepeak;
}


								/* Calculate coefficients for interpolation filter - lowpass filter after upsampling */
double *
calc_lp_os_coeffs (int samplerate, int os_factor, int taps)
{

  double *coeff_vector = calloc (taps, sizeof (double));	// will be freed from function that free TruePeak Context

  double nearzero = 0.000001;

  for (int i = 0; i < taps; i++)
    {
      double m = i - (taps - 1) / 2.0;
      double c = 1.0;
      if (fabs (m) >= nearzero)
	{
	  c = sin (m * M_PI / os_factor) / (m * M_PI / os_factor);
	  /* Should windowing be applied? Yes!! */
	  c = c * (0.5 * (1 - cos (2 * M_PI * i / (taps - 1))));
	  // To populate a polyphase I should do that otherwise
	  if (fabs (c) >= nearzero)
	    {
	      coeff_vector[i] = c;
	    }
	}

    }
  /*
     #ifdef DEBUG
     printf("Coefficients for low pass filter with %d taps\n", taps);
     for (int i = 0; i < taps; i++) {
     printf("%.8f\n", coeff_vector[i]);
     }
     #endif
   */
  return coeff_vector;
}



/*

   Things to do:
   - rework autotools and autoconfig to be sure it works on multiple platforms: see also problem with usability presence.
   - wouldn't it be better to use cmake  instead of autotools?
   - also depending on ffmpeg version the files should be opened differently
   - Look at the question of the first few milliseconds of DI output
   - last change because of ffmpeg 4.0. Also test actual software on older ffmpeg installations.
   - Look at automatic gate switching and threshold adjustments for leqm DI
   - ffmpeg and sndfile give differents results for DI. Look at that
   - add other Dialoque Intelligence Variants?
   - Look at the notes from Zurich
   - implement overlapping by percentual also for DolbyDI integration (render it modifiable)
   - implement level gating indipendent from dolby DI (at present can only be combined with dolbydi)
   - implement the A and C weightings for Leq. I could even switch from one or the other and various combinations
   - I should tend at using a single worker_function ( apart from di_worker_function for preprocessing )
   - delete workingbranch and the other no more used on gitlab and github
   - train tensorflow model for DI (I could use trailer for the training, automate data extraction)
   x add switch for printing out detailed DI information
   - report discarded samples for LKFS measurement
   -/ implement TruePeak according to ITU - POC is done, much in need of speed up through polyphase filters
   x implement Leq(M) according to LKFS logic with and without gating - worker_function_gated2
   - test audio data with only a single subdivision (single full gating Block at the end instead of say four)
   - LKFS premultiply a and b coefficients of the two filters instead of applying sequencially the 2 filters - in principle is clear
   x Found out that with normal default measurement the program measure slightly differently with changing number of cpus: this is probably because of deserialization. This must be cured
   x Also logging seems to be dependent on number of cpus: for ex. if 2 cpus buffersize interval ist doubled, with 3 tripled and so on.
   x Buffer B in LKFS must be scorporated from worker functions for problematic signalling of serialization and no present speed up with multiple cpus.
   
STRANGE THINGS HAPPENING:

- opening files with SNDFILE or FFMPEG produces the same LKFS and Leq(M) but different DI measurement including speech percentage !! Why?
Well, as DI has an automatic gaining system it is quite resilient with regard to changes of gain but very sensitive to the smallest changes 
in values in lower significant digits. As FFMPEG and SNDFILE do not return exactly the same normalized values (differences after the 9th oder 10th 
decimal digit) Dolby DI produces different results if opening files with sndfile or ffmpeg. So if this is the complete explication, implementing
a scaling independent from the used library should solve the issue.

- Leq(M), LKFS and also LKFS(DI) measurements scales perfectly (Leq(M) and LKFS) or almost (LKFS(DI), but Leq(M,LG) seems to half-scale, meaning that a reduction in 6 dB produces a reduction of only 3 dB. Why?
This should be an error somewhere.

- LKFS(DI) does not seems to pass the conformance test, as it measure 25.12 instead of 24 with ffmpeg and 25.21 with sndfile. But I doubt that this is possibly true because 23.1 is LKFS and there are passages without dialogue in the conformance test file, so how can this possibly be true? It could be that the problem has also to do with the metode used for opening wav files and normalizing (this could be checked to track down also the different dialogue percentages with different audio libraries.
*/
