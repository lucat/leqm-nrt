/*
    leqm-nrt is a  non-real-time implementation 
    of Leq(M) measurement according to ISO 21727:2004(E)
    "Cinematography -- Method of measurement of perceived
    loudness of motion-picture audio material"

    Copyright (C) 2011-2013, 2017-2019 Luca Trisciani

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

#define FFMPEG
//#define SNDFILELIB
#ifdef FFMPEG
#include <stdio.h>
#include <libavutil/avassert.h>
#include <libavcodec/avcodec.h>
#include <libavformat/avformat.h>
#include <libavutil/avutil.h>
#elif defined SNDFILELIB
#include <sndfile.h>
#include <samplerate.h>
#endif

// Version 0.0.20 (C) Luca Trisciani 2011-2013, 2017-2019
// Tool from the DCP-Werkstatt Software Bundle



// COMPILATION
// compile for DEBUG with gcc -g -DEBUG -I/usr/include/ffmpeg -lfftw3 -lm -lpthread -lrt  -lavformat -lavcodec -lavutil -o leqm-nrt leqm-nrt.c

//otherwise  gcc -I/usr/include/ffmpeg -lm -lpthread -lrt -lavformat -lavcodec -lavutil -o leqm-nrt leqm-nrt.c


//this to do:
/*

- with the github version and SNDFILELIB I get no memory leaks, so if I strip down of FFMPEG I should get the same here


 */



//#define DEBUG

#define VERSION 19

#ifdef SNDFILELIB
SRC_DATA src_data;
#endif

struct Sum {
  double csum; // convolved sum
  double sum; // flat sum
    int nsamples;
  double cmean; //convolved mean
    double mean;
    double leqm;
  double rms;
};

struct WorkerArgs {
  double * argbuffer;
  int nsamples;
  int nch;
  int npoints;
  double * ir;
  struct Sum * ptrtotsum;
  double * chconf;
  int shorttermindex;
  double * shorttermarray;
  int leqm10flag;
  #ifdef SNDFILELIB
  double src_output;
  #endif
};




int equalinterval( double * freqsamples, double * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints);
int equalinterval2( double freqsamples[], double * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile);
int convloglin(double * in, double * out, int points);
double convlinlog_single(double in);
double convloglin_single(double in);
int convolv_buff(double * sigin, double * sigout, double * impresp, int sigin_dim, int impresp_dim);
double inputcalib (double dbdiffch);
int rectify(double * squared, double * inputsamples, int nsamples);
int accumulatech(double * chaccumulator, double * inputchannel, int nsamples);
int sumsamples(struct Sum * ts, double * inputsamples, double * cinputsamples, int nsamples);
int meanoverduration(struct Sum * oldsum);
void  inversefft1(double * eqfreqresp, double * ir, int npoints);
void  inversefft2(double * eqfreqresp, double * ir, int npoints);
void * worker_function(void * argfunc);
void logleqm(FILE * filehandle, double featuretimesec, struct Sum * oldsum);
double sumandshorttermavrg(double * channelaccumulator, int nsamples);
double logleqm10(FILE * filehandle, double featuretimesec, double longaverage);
#ifdef FFMPEG
int transfer_decoded_data(AVFrame * ptr_frame, struct WorkerArgs ** dptrWorkerArgs, int w_id,  AVCodecContext* codecCon);
int transfer_decoded_samples(AVFrame * ptr_frame, double * buf, AVCodecContext* codecCon, int buffersizs, int * doublesample_index);
int transfer_remaining_decoded_samples(AVFrame * ptr_frame, double *bufremain, AVCodecContext* codecCon, int nxtsmpl, int * doublesample_index);
#endif
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;



#ifdef FFMPEG


static double get(uint8_t *a[], int ch, int index, int ch_count, enum AVSampleFormat f){
    const uint8_t *p;
    if(av_sample_fmt_is_planar(f)){
        f= av_get_alt_sample_fmt(f, 0);
        p= a[ch];
    }else{
        p= a[0];
        index= ch + index*ch_count;
    }

    switch(f){
    case AV_SAMPLE_FMT_U8 : return ((const uint8_t*)p)[index]/127.0-1.0;
    case AV_SAMPLE_FMT_S16: return ((const int16_t*)p)[index]/32767.0;
    case AV_SAMPLE_FMT_S32: return ((const int32_t*)p)[index]/2147483647.0;
    case AV_SAMPLE_FMT_FLT: return ((const float  *)p)[index];
    case AV_SAMPLE_FMT_DBL: return ((const double *)p)[index];
    default: av_assert0(0);
    }
}


void printAudioFrameInfo(const AVCodecContext* codecContext, const AVFrame* frame)
{
    // See the following to know what data type (unsigned char, short, float, etc) to use to access the audio data:
    // http://ffmpeg.org/doxygen/trunk/samplefmt_8h.html#af9a51ca15301871723577c730b5865c5
  printf("Audio frame info:\n");
  printf("  Sample count: %d\n", frame->nb_samples);
  printf("  Channel count: %d\n", codecContext->channels);
  printf("  Format: %s\n", av_get_sample_fmt_name(codecContext->sample_fmt));
  printf("  Bytes per sample: %d\n", av_get_bytes_per_sample(codecContext->sample_fmt));
  printf("  Is planar? %d\n", av_sample_fmt_is_planar(codecContext->sample_fmt));
  printf("frame->linesize[0] tells you the size (in bytes) of each plane\n");

  if (codecContext->channels > AV_NUM_DATA_POINTERS && av_sample_fmt_is_planar(codecContext->sample_fmt))
    {
      printf("The audio stream (and its frames) have too many channels to fit in\n" \
	     "frame->data. Therefore, to access the audio data, you need to use\n" \
             "frame->extended_data to access the audio data. It's planar, so\n" \
	     "each channel is in a different element. That is:\n" \
	     "  frame->extended_data[0] has the data for channel 1\n" \
	     "  frame->extended_data[1] has the data for channel 2\n" \
	     "  etc.\n");
    }
    else
    {
      printf("Either the audio data is not planar, or there is enough room in\n"\
	     "frame->data to store all the channels, so you can either use\n"\
	     "frame->data or frame->extended_data to access the audio data (they\n"\
	     "should just point to the same data).\n");
    }

	printf("If the frame is planar, each channel is in a different element.\n"\
	       "That is:\n"\
	       "  frame->data[0]/frame->extended_data[0] has the data for channel 1\n"\
	       "  frame->data[1]/frame->extended_data[1] has the data for channel 2\n"\
	       "  etc.\n");

       printf("If the frame is packed (not planar), then all the data is in\n"\
	      "frame->data[0]/frame->extended_data[0] (kind of like how some\n"\
              "image formats have RGB pixels packed together, rather than storing\n"\
              " the red, green, and blue channels separately in different arrays.\n");
}
#endif

int checkargvalue(const char * stringarg) {
  if (stringarg == NULL) {
      printf("Please provide required value after argument switch!\n");
      return 1;
    } else if (!(isdigit(stringarg[0]))) {
    if  ((strncmp(stringarg, "-", 1) == 0) && (isdigit(stringarg[1]))){
    return 0;
    } else {
      return 1;
    }
  } else {
      return 0;
}
}


int main(int argc, const char ** argv)
{
  int npoints = 64; // This value is low for precision. Calibration is done with 32768 point.
  int origpoints = 21; //number of points in the standard CCIR filter
  int samplingfreq; // this and the next is defined later taking it from sound file
  int bitdepth;
  // double normalizer;
  int timing = 0;
  struct timespec starttime;
  int fileopenstate = 0;
  int leqm10 = 0;
  int leqmlog = 0;
  #if defined __unix__ || defined  __APPLE__
  int numCPU = sysconf(_SC_NPROCESSORS_ONLN) - 1;
  #elif defined _WIN64 || defined _WIN32
  SYSTEM_INFO sysinfo;
  GetSystemInfo(&sysinfo);
  int numCPU = sysinfo.dwNumberOfProcessors - 1;
  #endif

  double * channelconfcalvector;
  channelconfcalvector = NULL;
  printf("leqm-nrt  Copyright (C) 2011-2013, 2017-2019 Luca Trisciani\nThis program comes with ABSOLUTELY NO WARRANTY; for details on command line parameters -help\nThis is free software, and you are welcome to redistribute it\nunder the GPL v3 licence.\nProgram will use 1 + %d slave threads.\n", numCPU);
  //SndfileHandle file;
#ifdef SNDFILELIB
  SNDFILE *file;
  file=NULL;
  SF_INFO sfinfo;
  memset(&sfinfo, 0, sizeof(sfinfo));
  int src_ratio = 8; //fix at present, could become a parameter
#elif defined FFMPEG
  av_register_all();
    AVFrame * frame = NULL;
    frame = av_frame_alloc();
    AVFormatContext* formatContext = NULL;
    AVCodec * cdc = NULL;
    AVStream* audioStream = NULL;
    AVCodecContext* codecContext = NULL;
    

#endif
    
  FILE *leqm10logfile;
  leqm10logfile = NULL;
  FILE *leqmlogfile;
  leqmlogfile = NULL;
  int buffersizems = 850; //ISO 21727:2004 do not contain any indication, TASA seems to indicate 1000, p. 8
  int buffersizesamples;
  	double tempchcal[128];
	int numcalread = 0;
	double longperiod = 10.0; //this is in minutes and correspond to the period for leqm10
	double threshold = 80.0; //this is the threshold for the Allen metric
	double * shorttermaveragedarray;
	shorttermaveragedarray = NULL;
	int numbershortperiods = 0;
	#ifdef FFMPEG
	int realnumbershortperiods = 0;
	#endif
	int parameterstate = 0;
	int leqnw = 0;

	char soundfilename[2048];
	// This is a requirement of sndfile library, do not forget it.

	const char helptext[] = "Order of parameters after audio file is free.\nPossible parameters are:\n--convpoints <integer number> \tNumber of interpolation points for the filter.\n\t\t\t\tDefault 64.\n--numcpus <integer number> \tNumber of slave threads to speed up operation.\n--timing \t\t\tFor benchmarking speed.\n--chconfcal <db correction> <db correction> <etc. so many times as channels>\n--logleqm10\t\t\t(will also print Allen metric as output)\n--threshold <leqm>\t\tThreshold used for Allen metric (default 80)\n--longperiod <minutes>\t\tLong period for leqm10 (default 10)\n--logleqm\n--buffersize <milliseconds>\nUsing:\ngnuplot -e \"plot \\\"logfile.txt\\\" u 1:2; pause -1\"\nit is possible to directly plot the logged data.\n";

	
  if (argc == 1)
    { 
      printf(helptext);
      printf("Please indicate a sound file to be processed.\n");
      return 0;
  }

    
    for (int in = 1; in < argc;) {
      if ((!(strncmp(argv[in], "-", 1) == 0)) && (argv[in] != NULL)) {

	  #ifdef SNDFILELIB
		if (fileopenstate == 0) {
	  if(! (file = sf_open(argv[in], SFM_READ, &sfinfo))) {
	    printf("Error while opening audio file, could not open  %s\n.", argv[in]);
	    puts(sf_strerror(NULL));
	    return 1;
	  }


	  strcpy(soundfilename, argv[in]);
	     fileopenstate = 1;
	     printf("Opened file: %s\n", argv[in]);
	     printf("Sample rate: %d\n", sfinfo.samplerate);
	     printf("Channels: %d\n", sfinfo.channels);
	     printf("Format: %d\n", sfinfo.format);
	     printf("Frames: %d\n", (int) sfinfo.frames);
	     channelconfcalvector = malloc(sizeof(double) * sfinfo.channels);
	     in++;
	     continue;
	} else {
	  free(channelconfcalvector);
	  channelconfcalvector = NULL;
	  return 0;
	}
      

	  #elif defined FFMPEG
			if (fileopenstate == 0) {
   if (avformat_open_input(&formatContext, argv[in], NULL, NULL) != 0)
    {
      //av_free(frame);
      av_frame_free(&frame);
        printf("Error opening the file\n");
        return 1;
    }
    //fileopenstate = 1;
    
    if (avformat_find_stream_info(formatContext, NULL) < 0)
    {
      //av_free(frame);
            av_frame_free(&frame);
        avformat_close_input(&formatContext);
        printf("Error finding the stream info\n");
        return 1;
    }

    // Find the audio stream
    //AVCodec* cdc = NULL;
    int streamIndex = av_find_best_stream(formatContext, AVMEDIA_TYPE_AUDIO, -1, -1, &cdc, 0);
    if (streamIndex < 0)
    {
      //av_free(frame);
            av_frame_free(&frame);
        avformat_close_input(&formatContext);
        printf("Could not find any audio stream in the file\n");
        return 1;
    }

    audioStream = formatContext->streams[streamIndex];
    codecContext = audioStream->codec;
    codecContext->codec = cdc;

    if (avcodec_open2(codecContext, codecContext->codec, NULL) != 0)
    {
      //      av_free(frame);
            av_frame_free(&frame);
        avformat_close_input(&formatContext);
        printf("Couldn't open the context with the decoder\n");
        return 1;
    }

	  strcpy(soundfilename, argv[in]);
	      fileopenstate = 1;
	    printf("This stream has %d ", codecContext->channels);
    printf(" channels and a sample rate of %d ", codecContext->sample_rate);
    printf(" Hz\n"); 
    printf("The data is in the format %s\n", av_get_sample_fmt_name(codecContext->sample_fmt)); 

	     channelconfcalvector = malloc(sizeof(double) *codecContext->channels);
	     in++;
	     continue;
			} else { /*if (fileopenstate == 0) */
	  free(channelconfcalvector);
	  channelconfcalvector = NULL;
	  printf("Please specify and input audio file.\n");
	  return 0;
	}

          
#endif
			continue;
      }   /*(!(strncmp(argv[in], "-", 1) == 0 */


      if (strcmp(argv[in], "--chconfcal") == 0) {
	/* as the order of parameter is free I have to postpone 
	   the check for consistency with the number of channels.
	   So first create a temporary array, whose number of element will be checked after 
	   the parsing of the command line parameters is finished. 
	   The calibration will be expressed in dB on the command line and converted to multiplier 
	   here so that it can be stored as a factor in the channelconfcalvector.
	*/

	in++;
	for (;;)  {
	if (in < argc) {
	  //if (!(strncmp(argv[in], "-", 1) == 0)) { //changed this to allow negative numbers
	    if (!(strncmp(argv[in], "-", 1) == 0) || isdigit(argv[in][1])) {
	  tempchcal[numcalread++]=atof(argv[in++]);
	  } else break;
	} else break;
	
	} //for
	continue;
      }
 
      if (strcmp(argv[in], "--convpoints") == 0)  {
       if (checkargvalue(argv[in + 1])) return 1;
	     npoints = atoi(argv[in + 1]);
	     in+=2;
	     printf("Convolution points sets to %d.\n", npoints);
	     continue;
	
      }
	      if (strcmp(argv[in], "--numcpus") == 0) {
		if (checkargvalue(argv[in + 1])) return 1;
		numCPU= atoi(argv[in + 1]);
	     in+=2;
	     printf("Number of threads manually set to %d. Default is number of cores in the system minus one.\n", numCPU);
	     continue;
	
      }
	      if (strcmp(argv[in], "--timing") == 0) {
		timing = 1;
	     in++;
	     printf("Execution time will be measured.\n");
	     continue;
	
      }
	            if (strcmp(argv[in], "--version") == 0) {
	     in++;
	     printf("This is leqm-nrt version 0.%d.\n", VERSION);
	     return 0;
	
      }
		    if (strcmp(argv[in], "--help") == 0) {
	     in++;
	     printf(helptext);
	     return 0;
	
      }

	      	      if (strcmp(argv[in], "--logleqm10") == 0) {
		leqm10 = 1;
	     in++;
	     printf("Leq(M)10 data will be logged to the file leqm10.txt\n");
	     continue;
	
      }
		      	      	      if (strcmp(argv[in], "--logleqm") == 0) {
		leqmlog = 1;
	     in++;
	     printf("Leq(M) data will be logged to the file leqmlog.txt\n");
	     continue;
	
      }

	     	      	      	      if (strcmp(argv[in], "--leqnw") == 0) {
		leqnw = 1;
	     in++;
	     printf("Leq(nW) - unweighted -  will be outputted.\n");
	     continue;
	
      }

				        if (strcmp(argv[in], "--buffersize") == 0) {
					  if (checkargvalue(argv[in + 1])) return 1;
		buffersizems = atoi(argv[in + 1]);
	     in+=2;
	     printf("Buffersize will be set to %d milliseconds.\n", buffersizems);
	     continue;
	
      }
					if (strcmp(argv[in], "--threshold") == 0) {
					  if (checkargvalue(argv[in + 1])) return 1;
		threshold = atof(argv[in + 1]);
	     in+=2;
	     printf("Threshold for Allen metric set to %f Leq(M).\n", threshold);
	     continue;
	
      }
					if (strcmp(argv[in], "--longperiod") == 0) {
					   if (checkargvalue(argv[in + 1])) return 1;
		longperiod = atof(argv[in + 1]);
	     in+=2;
	     printf("Longperiod for Leq(M)X set to %f minutes.\n", longperiod);
	     continue;
	
      }

					if (parameterstate==0) {
					  printf("The command line switch you typed is not valid: %s\n", argv[in]);
					  break;
					}

      
    } /* for (int in=1; in < argc;) */
// Open audio file

//postprocessing parameters
#ifdef SNDFILELIB
    if (numcalread == sfinfo.channels) {
      for (int cind = 0; cind < sfinfo.channels; cind++) {
	channelconfcalvector[cind] = convloglin_single(tempchcal[cind]);
	
      }
    }
    else if ((numcalread == 0) && (sfinfo.channels == 6)) {
	double conf51[6] = {0, 0, 0, 0, -3, -3};
	for (int cind = 0; cind < sfinfo.channels; cind++) {
	  channelconfcalvector[cind] = convloglin_single(conf51[cind]);
	}
	printf("Using input channel calibration for 5.1 configuration:\n0 0 0 0 -3 -3\n");
    }
        else if ((numcalread == 0) && (sfinfo.channels == 8)) {
	  double conf71[8] = {0, 0, 0, 0, -3, -3, -3, -3};
	for (int cind = 0; cind < sfinfo.channels; cind++) {
	  channelconfcalvector[cind] = convloglin_single(conf71[cind]);
	}
	printf("Using input channel calibration for 7.1 configuration:\n0 0 0 0 -3 -3 -3 -3\n");

    }
#elif defined FFMPEG

    if (numcalread == codecContext->channels) {
      for (int cind = 0; cind < codecContext->channels; cind++) {
	channelconfcalvector[cind] = convloglin_single(tempchcal[cind]);
      }
    }
    else if ((numcalread == 0) && (codecContext->channels == 6)) {
	double conf51[6] = {0, 0, 0, 0, -3, -3};
	for (int cind = 0; cind < codecContext->channels; cind++) {
	  channelconfcalvector[cind] = convloglin_single(conf51[cind]);
	}
	printf("Using input channel calibration for 5.1 configuration:\n0 0 0 0 -3 -3\n");	
    }
    else if ((numcalread == 0) && (codecContext->channels == 8)) {
      double conf71[8] = {0, 0, 0, 0, -3, -3, -3, -3};
	for (int cind = 0; cind < codecContext->channels; cind++) {
	  channelconfcalvector[cind] = convloglin_single(conf71[cind]);
	}
	printf("Using input channel calibration for 7.1 configuration:\n0 0 0 0 -3 -3 -3 -3\n");
    }
#endif
     else {

      printf("Either you specified a different number of calibration than number of channels in the file or you do not indicate any calibration and the program cannot infer one from the number of channels. Please specify a channel calibration on the command line.\n");

      free(channelconfcalvector);
      channelconfcalvector = NULL;
      return 0;
     }




    if (leqm10) {
      char tempstring[2048];
      strcpy(tempstring, soundfilename);
      strcat(tempstring, ".leqm10.txt");
      leqm10logfile = fopen(tempstring, "w");
      if (leqm10logfile == NULL) {
	printf("Could not open file to write log leqm10 data!\n");
      }
    }




    if (leqmlog) {
      char tempstring[2048];
      strcpy(tempstring, soundfilename);
      strcat(tempstring, ".leqmlog.txt");
      leqmlogfile = fopen(tempstring, "w");
      if (leqmlogfile == NULL) {
	printf("Could not open file to write log leqm data!\n");	
      }
    }
    
      
  if (timing) {

  clock_gettime(CLOCK_MONOTONIC, &starttime);
  }
  



  
  // reading to a double or float buffer with sndfile take care of normalization
 /*
 static double  buffer[BUFFER_LEN]; // it seems this must be static. I don't know why
 */
 double * buffer;
 #ifdef FFMPEG
 double * remainbuffer;
 #endif
 // buffer = new double [BUFFER_LEN];
 //buffersizesamples = (sfinfo.samplerate*sfinfo.channels*buffersizems)/1000;

#ifdef SNDFILELIB
 if ((sfinfo.samplerate*buffersizems)%1000) {
#elif defined FFMPEG
   if ((codecContext->sample_rate*buffersizems)%1000) {
#endif
   printf("Please fine tune the buffersize according to the sample rate\n");
   //close file
   // free memory
   // write a function to do that 
   return 1;
 }


   #ifdef SNDFILELIB
  buffersizesamples = (sfinfo.samplerate*sfinfo.channels*buffersizems)/1000;
  buffer = malloc(sizeof(double)*buffersizesamples);

 samplingfreq = sfinfo.samplerate;

 
 if(leqm10) {
   
   //if duration < 10 mm exit

   double featdursec = sfinfo.frames / sfinfo.samplerate;
   if ((featdursec/60.0) < longperiod) {
     printf("The audio file is too short to measure Leq(m10).\n");
     return 0;
   }
   
 
   //how many short periods in overall duration //ATTENTION this calculation may return 0  for buffer sizes of less than 10 ms. Maybe cast to double!!
   int remainder = sfinfo.frames % (sfinfo.samplerate*buffersizems/1000);
   if (remainder == 0)  numbershortperiods = sfinfo.frames/(sfinfo.samplerate*buffersizems/1000); 
   else  numbershortperiods = sfinfo.frames/(sfinfo.samplerate*buffersizems/1000) + 1;
  
   //allocate array
   shorttermaveragedarray = malloc(sizeof(*shorttermaveragedarray)*numbershortperiods);
 } /* if (leqm10) */

#elif defined FFMPEG
  buffersizesamples = (codecContext->sample_rate*codecContext->channels*buffersizems)/1000;
  buffer = malloc(sizeof(double)*buffersizesamples);
  remainbuffer = malloc(sizeof(double)*buffersizesamples);
 samplingfreq = codecContext->sample_rate;
 // postpone this because I cannot get duration or number of frames
 // but I cannot postpone this!!
 //it seems I cannot get total number of audio frames in ffmpeg so I will simply allocate enough
 //memory for a 5 hours feature, ok?
 if (leqm10) {
   /* //This check is not possible here with ffmpeg
   double featdursec = formatContext->duration;
   printf("Total duration is %d secs.", (int) formatContext->duration);
   if ((featdursec/60.0) < longperiod) {
     printf("The audio file is too short to measure Leq(m10).\n");
     return 0;
   }

   */

   //numbershortperiods = (int) (180000.00 / (((double) codecContext->sample_rate) * (double) buffersizems/1000.00) + 1); //this is wrong, because 180000 cannot be be number of frames, also why should we devide by the sample rate if 18000 is just seconds?
   numbershortperiods = (int) (18000.00 / ((double) buffersizems/1000.00) + 1);
   shorttermaveragedarray = malloc(sizeof(*shorttermaveragedarray)*numbershortperiods);
 }


 /* if (leqm10) */
 

 #endif

 //End opening audio file

  //ISO 21727:2004(E)
  // M Weighting
  double freqsamples[] = {31, 63, 100, 200, 400, 800, 1000, 2000, 3150, 4000, 5000, 6300, 7100, 8000, 9000, 10000, 12500, 14000, 16000, 20000, 31500};
  double freqresp_db[] = {-35.5, -29.5, -25.4, -19.4, -13.4, -7.5, -5.6, 0.0, 3.4, 4.9, 6.1, 6.6, 6.4, 5.8, 4.5, 2.5, -5.6, -10.9, -17.3, -27.8, -48.3};
  
  double * eqfreqresp_db;
  eqfreqresp_db = malloc(sizeof(*eqfreqresp_db)*npoints);

  double * eqfreqsamples;
  eqfreqsamples = malloc(sizeof(*eqfreqsamples)*npoints);
  double * eqfreqresp;
  eqfreqresp = malloc(sizeof(*eqfreqresp)*npoints);
  double * ir;
  ir = malloc(sizeof(*ir)*npoints*2);


// And what to do for floating point sample coding?
#ifdef SNDFILELIB
   switch(sfinfo.format & SF_FORMAT_SUBMASK) {
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
   printf("No known bitdepth! Exiting ...\n");
   return -1;
   }

#elif defined FFMPEG
   //sono qui
   bitdepth = av_get_exact_bits_per_sample(codecContext->codec_id);
   #ifdef DEBUG
   printf("ffmpeg report %d bitdepth for the file.", bitdepth);
   #endif

#endif

  
   equalinterval2(freqsamples, freqresp_db, eqfreqsamples, eqfreqresp_db, npoints, samplingfreq, origpoints, bitdepth);
  convloglin(eqfreqresp_db, eqfreqresp, npoints);

    #ifdef DEBUG
    for (int i=0; i < npoints; i++) {
      printf("%d\t%.2f\t%.2f\t%.2f\n", i, eqfreqsamples[i], eqfreqresp_db[i], eqfreqresp[i]);  
    }
    #endif
    
    inversefft2(eqfreqresp, ir, npoints);

// read through the entire file

   struct Sum * totsum;
    totsum = malloc(sizeof(struct Sum));
    totsum->csum = 0.0;
    totsum->sum = 0.0;
    totsum->nsamples = 0;
    totsum->cmean = 0.0;
    totsum->mean = 0.0; // Do I write anything here?
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
 // Main loop through audio file

 int worker_id = 0;
 pthread_t tid[numCPU];
 struct WorkerArgs ** WorkerArgsArray;
 WorkerArgsArray = malloc(sizeof(struct WorkerArgs *)*numCPU);
 int staindex = 0; //shorttermarrayindex

#ifdef SNDFILELIB

 src_data.end_of_input = 0;/* Set this later. */

	/* Start with zero to force load in while loop. */
	src_data.input_frames = 0 ;
	src_data.data_in = buffer ;

	src_data.src_ratio = src_ratio ; //this could be done a separate parameter

	src_data.data_out = src_output ;
	src_data.output_frames = BUFFER_LEN /channels ;

 
 while((samples_read = sf_read_double(file, buffer, buffersizesamples)) > 0) {
   
#elif defined FFMPEG
  
    AVPacket readingPacket;
    av_init_packet(&readingPacket);
    int data_size = 0;
    int copiedsamples = 0; //also pointer to position  wherein to copy into the buffer

      while (av_read_frame(formatContext, &readingPacket) == 0)
    {

    


           if (readingPacket.stream_index == audioStream->index)
        {
            AVPacket decodingPacket = readingPacket;

            // Audio packets can have multiple audio frames in a single packet
            while (decodingPacket.size > 0)
            {
                // Try to decode the packet into a frame
                // Some frames rely on multiple packets, so we have to make sure the frame is finished before
                // we can use it
                int gotFrame = 0;
                int result = avcodec_decode_audio4(codecContext, frame, &gotFrame, &decodingPacket);

                if (result >= 0 && gotFrame)
                {
                    decodingPacket.size -= result;
                    decodingPacket.data += result;

                    // We now have a fully decoded audio frame
		    #ifdef DEBUG
                    printAudioFrameInfo(codecContext, frame);
		    #endif
		    
		    //here goes the copying to the multithreaded buffer
		    //but this buffer is in milliseconds and could be less oder more
		    //probably greater than the single frame

		    //copy as much samples as there is room in the buffer
		    
		      if (gotFrame)
		      {
			data_size = //this is in bytes
		      av_samples_get_buffer_size
		      (
		      NULL, 
		      codecContext->channels,
		      frame->nb_samples,
		      codecContext->sample_fmt,
		      1
		       );
 // check if copying frame data will overflow the buffer
// and copy only the samples necessary to completely fill the buffer
// save the remaining for the next copy
// write a function that uses get to change data to double and copy in worker buffer
		 
	  copiedsamples += (frame->nb_samples * frame->channels);
	  //         memcpy((char) ((void *) buffer), frame.data[0], data_size);	  
	  //transfer_decoded_data(frame, WorkerArgsArray, worker_id, codecContext);
	  // nextsample is next sample of frame to be copied
	  nextsample = transfer_decoded_samples(frame, buffer, codecContext, buffersizesamples, &dsindex);
	  #ifdef DEBUG
	  printf("Next sample index: %d\n", nextsample);
	  #endif
	  //// From here execute only if buffer is full of data

	  //if (copiedsamples >= buffersizesamples) {
	  while (copiedsamples >= buffersizesamples) {
	    realnumbershortperiods++; 
   WorkerArgsArray[worker_id]=malloc(sizeof(struct WorkerArgs));
   WorkerArgsArray[worker_id]->nsamples = buffersizesamples;
   WorkerArgsArray[worker_id]->nch = codecContext->channels;
   WorkerArgsArray[worker_id]->npoints=npoints;
   WorkerArgsArray[worker_id]->ir = ir;
   WorkerArgsArray[worker_id]->ptrtotsum = totsum;

   WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
   if (leqm10) {
   WorkerArgsArray[worker_id]->shorttermindex = staindex++;
   WorkerArgsArray[worker_id]->leqm10flag = 1;
   WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
   } else {
     WorkerArgsArray[worker_id]->shorttermindex = 0;
     WorkerArgsArray[worker_id]->leqm10flag = 0;
   }
   dsindex = 0;
	    // store rest in another buffer
	    if (copiedsamples > buffersizesamples) {
	    //copiedsamples = transfer_remaining_decoded_samples(frame, remainbuffer, codecContext, nextsample, &dsindex)*frame->channels;
	    	transfer_remaining_decoded_samples(frame, remainbuffer, codecContext, nextsample, &dsindex)*frame->channels;
	    } else {
	    	copiedsamples = 0;
	    }
	    //remaindertot = copiedsamples - buffersizesamples;
	    //copiedsamples = 0;
	  WorkerArgsArray[worker_id]->argbuffer = malloc(sizeof(double)*buffersizesamples); 
   memcpy(WorkerArgsArray[worker_id]->argbuffer, (void *)buffer, buffersizesamples*sizeof(double));
   if (copiedsamples > buffersizesamples) { //mistake here because I updated copiedsamples, it is no more the intended one
   //once buffer copied, copy rest if present in buffer for next cycle
     //but I have to update dsindex otherwise I will overwrite data on the next cycle
	 //ACHTUNG: the buffer must be the one of the next worker!
	   // ALSO: I have to account for n_samples if buffer is not full on the last buffer
   memcpy(buffer, remainbuffer, sizeof(double)*(copiedsamples - buffersizesamples)); //here I should copy only copiedsamples-buffersizesamples
   dsindex = copiedsamples  - buffersizesamples; //added
   copiedsamples = copiedsamples - buffersizesamples;
   //copiedsamples already initialized to the right value
   }
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_create(&tid[worker_id], &attr, worker_function, WorkerArgsArray[worker_id]);

   worker_id++;
   if (worker_id == numCPU) {
       worker_id = 0;
       //maybe here wait for all cores to output before going on
       for (int idxcpu = 0; idxcpu < numCPU; idxcpu++) {
       pthread_join(tid[idxcpu], NULL);
       free(WorkerArgsArray[idxcpu]->argbuffer);
       WorkerArgsArray[idxcpu]->argbuffer = NULL;
       free(WorkerArgsArray[idxcpu]);
       WorkerArgsArray[idxcpu] = NULL;
       }
              //simply log here your measurement it will be a multiple of your threads and your buffer
       if (leqmlog) {
	 meanoverduration(totsum); //update leq(m) until now and log it
       logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) codecContext->sample_rate), totsum );
	       } //endlog
   } //if (worker_id == numCPU)
   

	  } /// till here if buffer is full, but if not? It will never start doing calculations
		      } //if (gotFrame)
      if(data_size <= 0) {
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
    } //if result >= 0 && gotFrame
                else
                {
                    decodingPacket.size = 0;
                    decodingPacket.data = NULL;
                }



		
            } //while decodingPacket.size
        } // if readingPaket.stream

        // You *must* call av_free_packet() after each call to av_read_frame() or else you'll leak memory
        av_free_packet(&readingPacket);

    


    //end while worker_id
 /// End looping cores
    } /* while (av_read_frame(formatContext, &readingPacket) */ // main loop through file


      /* Adding insert for last samples at the end of the file - v. 28 */
      if (copiedsamples < buffersizesamples) {
	realnumbershortperiods++;
	worker_id = 0; //this should have been already like that
	WorkerArgsArray[worker_id]=malloc(sizeof(struct WorkerArgs));
	WorkerArgsArray[worker_id]->nsamples = copiedsamples; //is this correct?
   WorkerArgsArray[worker_id]->nch = codecContext->channels;
   WorkerArgsArray[worker_id]->npoints=npoints;
   WorkerArgsArray[worker_id]->ir = ir;
   WorkerArgsArray[worker_id]->ptrtotsum = totsum;

   WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
   if (leqm10) {
   WorkerArgsArray[worker_id]->shorttermindex = staindex++;
   WorkerArgsArray[worker_id]->leqm10flag = 1;
   WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
   } else {
     WorkerArgsArray[worker_id]->shorttermindex = 0;
     WorkerArgsArray[worker_id]->leqm10flag = 0;
   }
   dsindex = 0;
	    //remaindertot = copiedsamples - buffersizesamples;
	    //copiedsamples = 0;
	  WorkerArgsArray[worker_id]->argbuffer = malloc(sizeof(double)*copiedsamples); 
   memcpy(WorkerArgsArray[worker_id]->argbuffer, (void *)buffer, copiedsamples*sizeof(double));
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_create(&tid[worker_id], &attr, worker_function, WorkerArgsArray[worker_id]);

   
       pthread_join(tid[worker_id], NULL);
       free(WorkerArgsArray[worker_id]->argbuffer);
       WorkerArgsArray[worker_id]->argbuffer = NULL;
       free(WorkerArgsArray[worker_id]);
       WorkerArgsArray[worker_id] = NULL;
       if (leqmlog) {
	 meanoverduration(totsum); //update leq(m) until now and log it
       logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) codecContext->sample_rate), totsum );
	       } //endlog
   
      } /*if (copiedsamples < buffersizesamples) */




  
   // Some codecs will cause frames to be buffered up in the decoding process. If the CODEC_CAP_DELAY flag
    // is set, there can be buffered up frames that need to be flushed, so we'll do that
    if (codecContext->codec->capabilities & CODEC_CAP_DELAY)
    {
        av_init_packet(&readingPacket);
        // Decode all the remaining frames in the buffer, until the end is reached
        int gotFrame = 0;
        while (avcodec_decode_audio4(codecContext, frame, &gotFrame, &readingPacket) >= 0 && gotFrame)
        {
            // We now have a fully decoded audio frame
	  // so I also need to process this
	  #ifdef DEBUG
            printAudioFrameInfo(codecContext, frame);
	    #endif
        }
    }
   #endif
   #ifdef SNDFILELIB
   WorkerArgsArray[worker_id]=malloc(sizeof(struct WorkerArgs));
   WorkerArgsArray[worker_id]->nsamples = samples_read;
   WorkerArgsArray[worker_id]->nch = sfinfo.channels;
   WorkerArgsArray[worker_id]->npoints=npoints;
   WorkerArgsArray[worker_id]->ir = ir;
   WorkerArgsArray[worker_id]->ptrtotsum = totsum;

   WorkerArgsArray[worker_id]->chconf = channelconfcalvector;
   if (leqm10) {
   WorkerArgsArray[worker_id]->shorttermindex = staindex++;
   WorkerArgsArray[worker_id]->leqm10flag = 1;
   WorkerArgsArray[worker_id]->shorttermarray = shorttermaveragedarray;
   } else {
     WorkerArgsArray[worker_id]->shorttermindex = 0;
     WorkerArgsArray[worker_id]->leqm10flag = 0;
   }
   
   WorkerArgsArray[worker_id]->argbuffer = malloc(sizeof(double)*buffersizesamples);
   //   WorkerArgsArray[worder_id]->src_output = malloc(sizeof(double)*buffersizesamples); // this is for sample rate conversion, not yet used
   memcpy(WorkerArgsArray[worker_id]->argbuffer, buffer, samples_read*sizeof(double));
   

   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_create(&tid[worker_id], &attr, worker_function, WorkerArgsArray[worker_id]);

   worker_id++;
   if (worker_id == numCPU) {
       worker_id = 0;
       //maybe here wait for all cores to output before going on
       for (int idxcpu = 0; idxcpu < numCPU; idxcpu++) {
       pthread_join(tid[idxcpu], NULL);
       free(WorkerArgsArray[idxcpu]->argbuffer);
       WorkerArgsArray[idxcpu]->argbuffer = NULL;
       free(WorkerArgsArray[idxcpu]);
       WorkerArgsArray[idxcpu] = NULL;
       }
              //simply log here your measurement it will be a multiple of your threads and your buffer
       if (leqmlog) {
	 meanoverduration(totsum); //update leq(m) until now and log it
       logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) sfinfo.samplerate), totsum );
	       } //endlog
   }
   


    //end while worker_id
 /// End looping cores
 } // main loop through file


 
#endif
 
 //here I should wait for rest workers (< numcpu)
 //but I need to dispose of thread id.
 if (worker_id != 0) { // worker_id == 0 means the number of samples was divisible through the number of cpus
   for (int idxcpu = 0; idxcpu < worker_id; idxcpu++) { //worker_id is at this point one unit more than threads launched
     pthread_join(tid[idxcpu], NULL);
          free(WorkerArgsArray[idxcpu]->argbuffer);
     WorkerArgsArray[idxcpu]->argbuffer = NULL;
     free(WorkerArgsArray[idxcpu]);
     WorkerArgsArray[idxcpu] = NULL;
   }
        //also log here for a last value
       if (leqmlog) {
	 meanoverduration(totsum); //update leq(m) until now and log it
	 #ifdef SNDFILELIB
       logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) sfinfo.samplerate), totsum );
       #elif defined FFMPEG
      logleqm(leqmlogfile, ((double) totsum->nsamples)/((double) codecContext->sample_rate), totsum );
       #endif
	       } //endlog  
 }
 // mean of scalar sum over duration
 
 meanoverduration(totsum);
 if (leqnw) {
 printf("Leq(noW): %.4f\n", totsum->rms); // Leq(no Weighting)
 }
 printf("Leq(M): %.4f\n", totsum->leqm);

  if(timing) {
   struct timespec stoptime;
   long stoptimenanoseconds;
   long executionnanoseconds;
   clock_gettime(CLOCK_MONOTONIC, &stoptime);
   
   if (stoptime.tv_nsec < starttime.tv_nsec) {
     stoptimenanoseconds = 1000000000 + stoptime.tv_nsec;
   } else {
     stoptimenanoseconds = stoptime.tv_nsec;
   }
   executionnanoseconds = stoptimenanoseconds - starttime.tv_nsec;
   printf("Total execution time is %.6f seconds\n", ((double) stoptime.tv_sec) - ((double) starttime.tv_sec) + ((double) executionnanoseconds / 1000000000.00));
 }


 if (leqm10) {
   //count shorttimeperiods through iterations!
   printf("Number short period buffers is: %d.\n", realnumbershortperiods);
   double duration = ((double) realnumbershortperiods) * ((double) buffersizems / 1000.0);
 //add remainder in samples!
 //this is not precise, because the last buffer will not be full
   //Take the array with the short term accumulators
   //double interval = 10.0;
   //create a rolling average according to rolling interval
   int rollint; // in short 10*60 = 600 sec 600/0.850 

   //how many element of the array to consider for the rollint?
   //that is how many buffersizems in the interval - interval could be parameterized(?)
   double tempint = 60.0 * longperiod / (((double) buffersizems) /1000.0); 
   rollint = (int) tempint;
   //dispose of the rest
   if (tempint - ((double) rollint) > 0) {
     rollint += 1;
   }

  
   printf("Total duration in minutes is %.0f.\n", duration/60.0);
   if ((duration/60.0) < longperiod) {
     printf("The audio file is too short to measure Leq(m10).\n");
   fclose(leqm10logfile);
   free(shorttermaveragedarray);
   //free(allenmetricarray);
   shorttermaveragedarray = NULL;
   //allenmetricarray = NULL;
     return 0; //but if I really want to exit here I should free memory
   }

  


   //numbershortperiods = (int) (180000.00 / (((double) codecContext->sample_rate) * (double) buffersizems/1000.00) + 1); //this is wrong, because 180000 cannot be be number of frames, also why should we devide by the sample rate if 18000 is just seconds?
   
   //shorttermaveragedarray = malloc(sizeof(*shorttermaveragedarray)*numbershortperiods);
 

   
   double * allenmetricarray = malloc (sizeof(*allenmetricarray)*(realnumbershortperiods -rollint));

   
   //two loops
   //external loop
   int indexlong = 0;
   double temp_leqm10 = 0.0;
   while(indexlong < (realnumbershortperiods - rollint)) {

     double accumulator = 0;
     //internal loop
     double averagedaccumulator = 0;
     for (int indexshort = 0; indexshort < rollint; indexshort++) {
       
       accumulator += shorttermaveragedarray[indexshort+indexlong];
     } //end internal loop
     averagedaccumulator = accumulator/((double) rollint);
     temp_leqm10 = logleqm10(leqm10logfile, ((double) (indexlong+rollint)) * ((double) buffersizems / 1000.0), averagedaccumulator);
     if (temp_leqm10 > threshold) { //See Allen article, this seems quite high...
       allenmetricarray[indexlong]=temp_leqm10; 
     } else {
       allenmetricarray[indexlong] = 0.0;
     }
     indexlong++;
   } //end external loop

   double thresholdedsum = 0;
   for (int i = 0; i < (realnumbershortperiods - rollint); i++) {
     thresholdedsum += allenmetricarray[i];
   }
   //printf("Allen Metric: %d", (int) (thresholdedsum / ((double) numbershortperiods)); // But Ioan Allen seems to require minutes as unites.
   printf("Allen metric: %d.\n", (int) (thresholdedsum / (duration/60.0))); // But Ioan Allen seems to require minutes as unites. But considering that the buffers are set to 750 ms it will be essentially the same, simply spreaded out times 80.
   fclose(leqm10logfile);
   free(shorttermaveragedarray);
   free(allenmetricarray);
   shorttermaveragedarray = NULL;
   allenmetricarray = NULL;
 }

 
 if (leqmlog) {

   fclose(leqmlogfile);
 }
 
 free(eqfreqsamples);
 eqfreqsamples = NULL;
  free(eqfreqresp_db);
 eqfreqresp_db=NULL;
 free(eqfreqresp);
 eqfreqresp = NULL;
 free(ir);
 ir = NULL;
 free(channelconfcalvector);
 channelconfcalvector = NULL;
 free(WorkerArgsArray);
 WorkerArgsArray = NULL;
 
 free(totsum);
 totsum = NULL;
 free(buffer);
 buffer=NULL;
 #ifdef FFMPEG
 free(remainbuffer);
 remainbuffer=NULL;
 #endif
 #ifdef SNDFILELIB
     
 sf_close(file);
 #elif defined FFMPEG
   // Clean up!
    //av_free(frame);
          av_frame_free(&frame);
    avcodec_close(codecContext);
    avformat_close_input(&formatContext);
 #endif
   return 0;
 }
       


    

void * worker_function(void * argstruct) {

  struct WorkerArgs * thisWorkerArgs = (struct WorkerArgs *) argstruct;

   double * sumandsquarebuffer;
   double * csumandsquarebuffer;
  double * chsumaccumulator_norm;
  double * chsumaccumulator_conv;
  

  sumandsquarebuffer = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));
  

  csumandsquarebuffer = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_norm = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));

  chsumaccumulator_conv = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));



  for (int i = 0; i < thisWorkerArgs->nsamples / thisWorkerArgs->nch; i++) {
    sumandsquarebuffer[i] = 0.0;
    csumandsquarebuffer[i] = 0.0;
  chsumaccumulator_norm[i] = 0.0;
  chsumaccumulator_conv[i] = 0.0;
  }


  
  for (int ch = 0; ch < thisWorkerArgs->nch; ch++) {

    double * normalizedbuffer;
    double * convolvedbuffer;


    normalizedbuffer = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));

    convolvedbuffer = malloc(sizeof(double)*(thisWorkerArgs->nsamples / thisWorkerArgs->nch));
    

    for (int n=ch, m= 0; n < thisWorkerArgs->nsamples; n += thisWorkerArgs->nch, m++) {
     // use this for calibration depending on channel config for ex. chconf[6] = {1.0, 1.0, 1.0, 1.0, 0.707945784, 0.707945784} could be the default for 5.1 soundtracks
      //so not normalized but calibrated
   normalizedbuffer[m] = thisWorkerArgs->argbuffer[n]*thisWorkerArgs->chconf[ch]; //this scale amplitude according to specified calibration

   
 }

 //convolution
 convolv_buff(normalizedbuffer, convolvedbuffer, thisWorkerArgs->ir, thisWorkerArgs->nsamples / thisWorkerArgs->nch, thisWorkerArgs->npoints * 2);
 //rectify, square und sum
 rectify(csumandsquarebuffer,convolvedbuffer, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
 rectify(sumandsquarebuffer,normalizedbuffer, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
 
 accumulatech(chsumaccumulator_norm, sumandsquarebuffer, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
 accumulatech(chsumaccumulator_conv, csumandsquarebuffer, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
 

 free(normalizedbuffer);
 normalizedbuffer= NULL;

 free(convolvedbuffer);
 convolvedbuffer=NULL;

 } // loop through channels

    //Create a function for this also a tag so that the worker know if he has to do this or not

  if (thisWorkerArgs->leqm10flag) {
    thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex] = sumandshorttermavrg(chsumaccumulator_conv, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
    #ifdef DEBUG
    printf("%d: %.6f\n", thisWorkerArgs->shorttermindex, thisWorkerArgs->shorttermarray[thisWorkerArgs->shorttermindex]);
    #endif
  }
  pthread_mutex_lock(&mutex);
  // this should be done under mutex conditions -> shared resources!
  sumsamples(thisWorkerArgs->ptrtotsum, chsumaccumulator_norm, chsumaccumulator_conv, thisWorkerArgs->nsamples / thisWorkerArgs->nch);
  pthread_mutex_unlock(&mutex);
  

  free(sumandsquarebuffer);
  sumandsquarebuffer=NULL;

  free(csumandsquarebuffer);
  csumandsquarebuffer=NULL;

  free(chsumaccumulator_norm);
  chsumaccumulator_norm=NULL;

  free(chsumaccumulator_conv);
  chsumaccumulator_conv=NULL;

  free(thisWorkerArgs->argbuffer);
  thisWorkerArgs->argbuffer = NULL;
  // the memory pointed to by this pointer is freed in main
  // it is the same memory for all worker
  // but it is necessary to set pointer to NULL otherwise free will not work later (really?)
  thisWorkerArgs->chconf = NULL;
 pthread_exit(0);

}


 #ifdef FFMPEG
 int transfer_decoded_data(AVFrame * ptr_frame, struct WorkerArgs ** dptrWorkerArgs, int w_id,  AVCodecContext* codecCon){
   int doublesample_index = 0; //this is to index the millisecond buffer of the worker
		for ( int ch_index = 0; ch_index < ptr_frame->channels; ch_index++) {

		  for (int smpl_index = 0; smpl_index< ptr_frame->nb_samples; smpl_index++) { //limit to buffer measure and return rest?
		      
		      dptrWorkerArgs[w_id]->argbuffer[doublesample_index++] = get(ptr_frame->data, ch_index, smpl_index, ptr_frame->channels, codecCon->sample_fmt);
		    // end print samples 
		    } //for nb_samples
		    } //for channels
		return ptr_frame->nb_samples;
 }

 // will return 0 or the index of next to be copied sample, if not enough room in buffer
 int transfer_decoded_samples(AVFrame * ptr_frame, double * buf, AVCodecContext* codecCon, int buffersizs, int * doublesample_index){
   // static int doublesample_index = 0; //this is to index the millisecond buffer of the worker
   //Yes the sample loop must be the external one, as I work with interleaved channels, not planar channels
   for (int smpl_index = 0; smpl_index< ptr_frame->nb_samples; smpl_index++) { //limit to buffer measure and return rest?
     for ( int ch_index = 0; ch_index < ptr_frame->channels; ch_index++) {
		    buf[(*doublesample_index)++] = get(ptr_frame->data, ch_index, smpl_index, ptr_frame->channels, codecCon->sample_fmt);
		      #ifdef DEBUG
		    printf("%0.5f\n", buf[(*doublesample_index)-1]);
		      #endif
		      if ((*doublesample_index) == buffersizs) {
			(*doublesample_index) = 0;
			//return (ptr_frame->nb_samples - smpl_index) * ptr_frame->channels;
			return smpl_index+1;
		      }
		    // end print samples 
		    } //for nb_samples
		    } //for channels
		return ptr_frame->nb_samples;
 }


 int transfer_remaining_decoded_samples(AVFrame * ptr_frame, double *bufremain, AVCodecContext* codecCon, int nxtsmpl, int * doublesample_index){
   //int doublesample_index = 0; //this is to index the millisecond buffer of the worker
   //Yes the sample loop must be the external one, as I work with interleaved channels, not planar channels
   for (int smpl_index = nxtsmpl; smpl_index< ptr_frame->nb_samples; smpl_index++) { //limit to buffer measure and return rest?
     for ( int ch_index = 0; ch_index < ptr_frame->channels; ch_index++) {
       bufremain[(*doublesample_index)++] = get(ptr_frame->data, ch_index, smpl_index, ptr_frame->channels, codecCon->sample_fmt);
		      #ifdef DEBUG
		    printf("%0.5f\n", bufremain[(*doublesample_index)-1]);
		      #endif		     
		    // end print samples 
		    } //for nb_samples
		    } //for channels
		return (ptr_frame->nb_samples) - nxtsmpl; //but what if also the second round at the same frame fill the buffer?
 }
 
#endif
 
//to get impulse response frequency response at equally spaced intervals is needed

int equalinterval( double * freqsamples, double  * freqresp, double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints) {
    double freq;
    // int findex = 0;
    // int rindex = 0;
    double pass = ((double) (samplingfreq >> 1)) / ((double) points);
    for (int ieq = 0, i = 0; ieq < points; ieq++) {
        freq = ieq*pass;
        eqfreqsamples[ieq] = freq;
	
        if ((freq == 0.0) || (freq < freqsamples[1])) { 
	  eqfreqresp[ieq] = freqresp[0];
            continue;
    } else {
        
        if ((freq >= freqsamples[i]) && (freq < freqsamples[i+1])) {
	  eqfreqresp[ieq] = ((freqresp[i+1] - freqresp[i])/(freqsamples[i+1] - freqsamples[i]))*(freq - freqsamples[i]) + freqresp[i];
        } else if (freq >=freqsamples[i+1]) {
            while(freq >= freqsamples[i+1]) {
                i++;
		if ((i + 1) >= origpoints) { 
		  break;
		}
            }
	    if ((i+1) < origpoints) {
            eqfreqresp[ieq] = ((freqresp[i+1] - freqresp[i])/(freqsamples[i+1] - freqsamples[i]))*(freq- freqsamples[i]) + freqresp[i];
	    } else {
	      eqfreqresp[ieq] = ((1 - freqresp[i])/(((double) (samplingfreq >> 1)) - freqsamples[i]))*(freq- freqsamples[i]) + freqresp[i];
	    }
        }
        }
    }
    return 0;
}





//the following is different from version 1 because interpolate between db and not linear. Conversion from db to lin must be done after.
//it is also different for the way it interpolates between DC and 31 Hz
// Pay attention that also arguments to the functions are changed
int equalinterval2( double freqsamples[], double  freqresp_db[], double * eqfreqsamples, double * eqfreqresp, int points, int samplingfreq, int origpoints, int bitdepthsoundfile) {
    double freq;


    //calculate miminum attenuation depending on the bitdeph (minus one), that is −6.020599913 dB per bit in eccess to sign
    double dcatt = ((double) (bitdepthsoundfile - 1))*(-6.020599913) + 20.00; //in dB
    //double dcatt = -90.3;
    double pass = ((double) (samplingfreq >> 1)) / ((double) points);
    for (int ieq = 0, i = 0; ieq < points; ieq++) {
        freq = ieq*pass;
        eqfreqsamples[ieq] = freq;
	
        if (freq == 0.0) {
	  eqfreqresp[ieq] = dcatt;
	} else if (freq < freqsamples[0]) { // this has a lot of influence on final Leq(M) value
	  eqfreqresp[ieq] = ((freqresp_db[0] - dcatt) / (freqsamples[0] - 0)) * freq + dcatt;
	  //eqfreqresp[ieq] = freqresp_db[0]; // Is this meaningful? Shouldn't I interpolate between 0 Hz and 31 Hz? Otherwise for DC I have -35.5 dB
            continue;
    } else {
        
        if ((freq >= freqsamples[i]) && (freq < freqsamples[i+1])) {
	  eqfreqresp[ieq] = ((freqresp_db[i+1] - freqresp_db[i])/(freqsamples[i+1] - freqsamples[i]))*(freq - freqsamples[i]) + freqresp_db[i];
        } else if (freq >=freqsamples[i+1]) {
            while(freq >= freqsamples[i+1]) {
                i++;
		if ((i + 1) >= origpoints) { 
		  break;
		}
            }
	    if ((i+1) < origpoints) {
            eqfreqresp[ieq] = ((freqresp_db[i+1] - freqresp_db[i])/(freqsamples[i+1] - freqsamples[i]))*(freq- freqsamples[i]) + freqresp_db[i];
	    } else {
	      eqfreqresp[ieq] = ((1 - freqresp_db[i])/(((double) (samplingfreq >> 1)) - freqsamples[i]))*(freq- freqsamples[i]) + freqresp_db[i];
	    }
        }
        }
    }
    return 0;
}






int convloglin(double * in, double * out, int points) {
  for (int i = 0; i < points; i++) {
    out[i] = powf(10, (in[i]/20.0));
  }

  return 0;
}


double convlinlog_single(double in) {
  double out;
    out = log(in)*20.0f;
  return out;
}


double convloglin_single(double in) {
  double out;
  out = powf(10, in/20.0f);
  return out;
}

// convolution

int convolv_buff(double * sigin, double * sigout, double * impresp, int sigin_dim, int impresp_dim) {


  double  sum = 0.0;
  for (int i = 0; i < sigin_dim; i++) {

    int m = i;
    for (int l = impresp_dim - 1; l >=0; l--,m++) {
      if (m >= sigin_dim) {
	m -= sigin_dim;
      }
      sum += sigin[m]*impresp[l];
    }
    sigout[i] = sum;
    sum=0.0;
    }
  return 0; 
  
}


void  inversefft2(double * eqfreqresp, double * ir, int npoints) {
  for (int n = 0; n < npoints; n++) {
    double parsum = 0.0;
    double partial = 0.0;
    
    for (int m = 1; m <= npoints -1; m++) {
      partial = cos(2.0*M_PI*((double) m)*( ( ((double) n) - ( ((double) npoints) * 2.0 -1 ) / 2 ) / ( ((double) npoints) * 2.0) ));
      parsum = parsum + eqfreqresp[m]*partial;
    }
    ir[n] = (eqfreqresp[0] + 2.0 * parsum)/((double) npoints * 2.0);
    #ifdef DEBUG
    printf("%.4f\n", ir[n]);
    #endif
  }
  for (int n = 0; n < npoints; n++) {
    ir[npoints+n] = ir[npoints-(n + 1)];
    #ifdef DEBUG
    printf("%.4f\n", ir[npoints+n]);
    #endif
  }
  
  
}

// scale input according to required calibration
// this could be different for certain digital cinema formats
double inputcalib (double dbdiffch) {
    
    double coeff = pow(10, dbdiffch/20);
    return coeff;
    
}

//rectify, square and sum
int rectify(double * squared, double * inputsamples, int nsamples){
  for (int i = 0; i < nsamples; i++) {
    squared[i] = (double) powf(inputsamples[i], 2);
    }
    return 0; 
    
}

int initbuffer(double * buffertoinit, int nsamples) {
  for (int i = 0; i < nsamples; i++) {
    buffertoinit[i] = 0.0;

  }
  return 0;
}

int accumulatech(double * chaccumulator, double * inputchannel, int nsamples) {
  for (int i = 0; i < nsamples; i++) {
    chaccumulator[i] += inputchannel[i];
  }
  return 0;
}

int sumsamples(struct Sum * ts, double * inputsamples, double * cinputsamples, int nsamples) {
  ts->nsamples += nsamples;
  for (int i=0; i < nsamples; i++) {
    ts->sum  += inputsamples[i];
    ts->csum += cinputsamples[i];
  }
  return 0;
  
}

int meanoverduration(struct Sum * oldsum) {
  oldsum->mean = pow(oldsum->sum / ((double) oldsum->nsamples), 0.500);
   oldsum->cmean = pow(oldsum->csum / ((double) oldsum->nsamples), 0.500);
   oldsum->rms = 20*log10(oldsum->mean) + 108.010299957;
   if (oldsum->rms < 0.0) {
     oldsum->rms = 0.0;
   }
   oldsum->leqm = 20*log10(oldsum->cmean) + 108.010299957;//
   if (oldsum->leqm < 0.0) {
     oldsum->leqm = 0.0;
   }

   /*
How the final offset is calculated without reference to a test tone:
P0 is the SPL reference 20 uPa

Reference SPL is RMS ! So 85 SPL over 20 uPa is 10^4.25 x 0.000020 = 0.355655882 Pa (RMS), 
but Peak value is 0.355655882 x sqr(2) = 0.502973372 that is 20 x log ( 0.502973372 / 0.000020) = 88.010299957

To that one has to add the 20 dB offset of the reference -20dBFS: 88.010299957 + 20.00 = 108.010299957 
   */
   /*But ISO 21727:2004(E) ask for a reference level "measured using an average responding meter". So reference level is not 0.707, but 0.637 = 2/pi
   */
return 0;
}

double sumandshorttermavrg(double * channelaccumulator, int nsamples) {
  double stsum = 0.0;
  for (int i=0; i < nsamples; i++) {
    stsum += channelaccumulator[i];
    
  }
  return stsum / (double) nsamples;
}

void logleqm(FILE * filehandle, double featuretimesec, struct Sum * oldsum) {

  fprintf(filehandle, "%.4f", featuretimesec);
  fprintf(filehandle, "\t");
  fprintf(filehandle, "%.4f\n", oldsum->leqm);
  

}

double logleqm10(FILE * filehandle, double featuretimesec, double longaverage) {
  double leqm10 = 20*log10(pow(longaverage, 0.500)) + 108.010299957;
  if (leqm10 < 0.0) {
    leqm10 = 0.0;
  }
  fprintf(filehandle, "%.4f", featuretimesec);
  fprintf(filehandle, "\t");
  fprintf(filehandle, "%.4f\n", leqm10);
  return leqm10;
}
