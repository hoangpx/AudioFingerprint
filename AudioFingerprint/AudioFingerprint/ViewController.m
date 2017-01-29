//
//  ViewController.m
//  AudioFingerprint
//
//  Created by Hoang Pham on 1/11/17.
//  Copyright © 2017 Hoang Pham. All rights reserved.
//

#import "ViewController.h"
#import <AudioToolbox/AudioToolbox.h>
#import "NSString+Util.h"
#import <CommonCrypto/CommonDigest.h>
#import <Accelerate/Accelerate.h>
#import <AVFoundation/AVFoundation.h>
#define CAP_SAMPLE_RATE 44100
#define CAP_CHANNELS 2
#define CAP_SAMPLE_SIZE sizeof(Float32) //4
AudioStreamBasicDescription const CAPAudioDescription = {
    .mSampleRate        = CAP_SAMPLE_RATE,
    .mFormatID          = kAudioFormatLinearPCM,
    .mFormatFlags       = kAudioFormatFlagIsSignedInteger,
    .mBytesPerPacket    = CAP_SAMPLE_SIZE * CAP_CHANNELS,
    .mFramesPerPacket   = 1,
    .mBytesPerFrame     = CAP_CHANNELS * CAP_SAMPLE_SIZE,
    .mChannelsPerFrame  = CAP_CHANNELS,
    .mBitsPerChannel    = 8 * CAP_SAMPLE_SIZE, //8 bits per byte
    .mReserved          = 0
};


@interface ViewController ()

@end

@implementation ViewController
CAPAudioOutput _audioOutput;
- (void)viewDidLoad {
    [super viewDidLoad];
    [self startProcessFile5];
//    [self startProcessFile2:&_audioOutput.player];
////    // Start audio
//    CAPStartAudioOutput(&_audioOutput);
}


- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

- (void)startProcessFile4 {
    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"Brad-Sucks--Total-Breakdown"
                                                              ofType:@"mp3"];
    NSData *data = [NSData dataWithContentsOfFile:audioFilePath];
    NSString *result = [self NSDataToHex:data];
    
}

static inline char itoh(int i) {
    if (i > 9) return 'A' + (i - 10);
    return '0' + i;
}

- (NSString *) NSDataToHex:(NSData *)data {
    NSUInteger i, len;
    unsigned char *buf, *bytes;
    
    len = data.length;
    bytes = (unsigned char*)data.bytes;
    buf = malloc(len*2);
    
    for (i=0; i<len; i++) {
        buf[i*2] = itoh((bytes[i] >> 4) & 0xF);
        buf[i*2+1] = itoh(bytes[i] & 0xF);
    }
    
    return [[NSString alloc] initWithBytesNoCopy:buf
                                          length:len*2
                                        encoding:NSASCIIStringEncoding
                                    freeWhenDone:YES];
}

- (void)startProcessFile5 {
    AudioStreamBasicDescription const audioDescription = {
        .mSampleRate        = 44100,
        .mFormatID          = kAudioFormatLinearPCM,
        .mFormatFlags       = kAudioFormatFlagIsSignedInteger,
        .mBytesPerPacket    = 8,
        .mFramesPerPacket   = 1,
        .mBytesPerFrame     = 8,
        .mChannelsPerFrame  = 2,
        .mBitsPerChannel    = 32,
        .mReserved          = 0
    };
    
    ExtAudioFileRef audioFile;
    OSStatus status = noErr;
    // Open file
    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"Brad-Sucks--Total-Breakdown"
                                                              ofType:@"mp3"];
//    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"sample"
//                                                                                         ofType:@"m4a"];
    NSURL *audioFileURL = [NSURL fileURLWithPath:audioFilePath];
    ExtAudioFileOpenURL((__bridge CFURLRef)audioFileURL, &audioFile);
    
    // Get files information
    AudioStreamBasicDescription fileAudioDescription;
    UInt32 size = sizeof(fileAudioDescription);
    ExtAudioFileGetProperty(audioFile,
                            kExtAudioFileProperty_FileDataFormat,
                            &size,
                            &fileAudioDescription);
    
    // Apply audio format
    ExtAudioFileSetProperty(audioFile,
                            kExtAudioFileProperty_ClientDataFormat,
                            sizeof(audioDescription),
                            &audioDescription);
    
    // Determine length in frames (in original file's sample rate)
    SInt64 fileLengthInFrames;
    size = sizeof(fileLengthInFrames);
    ExtAudioFileGetProperty(audioFile,
                            kExtAudioFileProperty_FileLengthFrames,
                            &size,
                            &fileLengthInFrames);

    // Calculate the true length in frames, given the original and target sample rates
    fileLengthInFrames = ceil(fileLengthInFrames * (audioDescription.mSampleRate / fileAudioDescription.mSampleRate));
    
    // Prepare AudioBufferList: Interleaved data uses just one buffer, non-interleaved has two
    int numberOfBuffers = audioDescription.mFormatFlags & kAudioFormatFlagIsNonInterleaved ? audioDescription.mChannelsPerFrame : 1;
    int channelsPerBuffer = audioDescription.mFormatFlags & kAudioFormatFlagIsNonInterleaved ? 1 : audioDescription.mChannelsPerFrame;
    int bytesPerBuffer = audioDescription.mBytesPerFrame * (int)fileLengthInFrames;
    
    AudioBufferList *bufferList = malloc(sizeof(AudioBufferList) + (numberOfBuffers-1)*sizeof(AudioBuffer));

    UInt32 readFrames = 0;
    float *audioTotal = malloc(fileLengthInFrames * sizeof(float));
    int index = 0;
    while (readFrames < fileLengthInFrames) {
        UInt32 framesLeftToRead = (UInt32)fileLengthInFrames - readFrames;
        UInt32 framesToRead = (4096 < framesLeftToRead) ? 4096 : framesLeftToRead;
        bufferList->mNumberBuffers = numberOfBuffers;
        bufferList->mBuffers[0].mData = calloc(bytesPerBuffer, 1);
        bufferList->mBuffers[0].mDataByteSize = bytesPerBuffer;
        bufferList->mBuffers[0].mNumberChannels = channelsPerBuffer;
        
        ExtAudioFileRead(audioFile, &framesToRead, bufferList);
        
        // store data
        SInt16 *inputFrames = (SInt16*)bufferList->mBuffers[0].mData;

        for(int i = 0; i < framesToRead; ++i) {
            audioTotal[index] = (float)inputFrames[i];
            index ++;
        }

        if ( framesToRead == 0 )
            break;
        
        readFrames += framesToRead;

    }
    
    ExtAudioFileDispose(audioFile);
//    free(bufferList->mBuffers[0].mData);
//    free(audioTotal);
//    free(bufferList);
    
    
    // process data
    int chunkSize = 512;
    int totalSize = index;
    int sampledChunkSize = totalSize/chunkSize;
    NSMutableArray *result = [[NSMutableArray alloc] init];
    NSMutableArray *highscores = [[NSMutableArray alloc] init];
    NSMutableArray *points = [[NSMutableArray alloc] init];
       
    for(int j = 0;j < sampledChunkSize; j++) {
        float *complexArray = malloc(sizeof(float) * chunkSize);
        NSMutableArray *tempResult = [[NSMutableArray alloc] init];
        for(int i = 0; i < chunkSize; i++) {
            complexArray[i] = audioTotal[(j*chunkSize)+i];
        }
        
        // Setup the length
        vDSP_Length log2n = log2f(chunkSize);
        
        // Calculate the weights array. This is a one-off operation.
        FFTSetup fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        
        // For an FFT, numSamples must be a power of 2, i.e. is always even
        int nOver2 = chunkSize/2;
        // I dont know purpose of this block
//         Populate *window with the values for a hamming window function
//        float *window = (float *)malloc(sizeof(float) * chunkSize);
//        vDSP_hamm_window(window, chunkSize, 0);
//        // Window the samples
//        vDSP_vmul(complexArray, 1, window, 1, complexArray, 1, chunkSize);
        
        // Define complex buffer
        COMPLEX_SPLIT A;
        A.realp = (float *) malloc(nOver2*sizeof(float));
        A.imagp = (float *) malloc(nOver2*sizeof(float));
        
        // Pack samples:
        // C(re) -> A[n], C(im) -> A[n+1]
        vDSP_ctoz((COMPLEX*)complexArray, 2, &A, 1, nOver2);
        
        //Perform a forward FFT using fftSetup and A
        //Results are returned in A
        vDSP_fft_zrip(fftSetup, &A, 1, log2n, FFT_FORWARD);
        
        //Convert COMPLEX_SPLIT A result to magnitudes
        float *amp = (float *) malloc(chunkSize*sizeof(float));
        amp[0] = A.realp[0]/(chunkSize*2);
        
        for(int i = 1; i < chunkSize; i++) {
            // get magnitudes
            amp[i] = sqrt(A.realp[i]*A.realp[i] + A.imagp[i]*A.imagp[i]);
            printf("%f ",amp[i]); // value from 0 to 617779.875000
            [tempResult addObject:[NSNumber numberWithFloat: amp[i]]];
        }
        NSMutableArray *temp = [[NSMutableArray alloc] init];
        NSMutableArray *temp1 = [[NSMutableArray alloc] init];
        for (int i =0; i < 5; i ++) {
            [temp addObject:[NSNumber numberWithFloat:0.0]];
            [temp1 addObject:[NSNumber numberWithFloat:0.0]];
        }
        [points addObject:temp1];
        [highscores addObject:temp];
        for (int freq = 40; freq < 300 ; freq++) {
            // Get the magnitude:
            float mag = ((NSNumber *)[tempResult objectAtIndex:freq]).floatValue;
            
            // Find out which range we are in:
            int index = [self getIndex:freq];
            
            // Save the highest magnitude and corresponding frequency:
            if (mag > ((NSNumber *)[[highscores objectAtIndex:j] objectAtIndex:index]).floatValue) {
                [[points objectAtIndex:j] replaceObjectAtIndex:index withObject: [NSNumber numberWithFloat:freq]];
                [[highscores objectAtIndex:j] replaceObjectAtIndex:index withObject:[NSNumber numberWithFloat:mag]];
            }
        }
        NSTimeInterval duration = chunkSize*j/44100;
//        printf("| time %f ", duration);
        long first = [[[points objectAtIndex:j] objectAtIndex:0] longValue];
        long second = [[[points objectAtIndex:j] objectAtIndex:1] longValue];
        long third = [[[points objectAtIndex:j] objectAtIndex:2] longValue];
        long forth = [[[points objectAtIndex:j] objectAtIndex:3] longValue];
        long h = [self hash:first p2:second p3:third p4:forth];
//        NSLog(@"%ld %ld %ld %ld ", first, second, third, forth);
        
        [result addObject:[NSNumber numberWithLong:h]];
        free(amp);
    }
    
    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString *documentsDirectory = [paths objectAtIndex:0];
    NSString *filePath = [NSString stringWithFormat:@"%@/sample.csv",documentsDirectory];
    [result writeToFile:filePath atomically:YES];

}

-(long) hash:(long) p1 p2:(long) p2 p3:(long) p3 p4:(long) p4 {
//    NSLog(@"%ld", p4);
    int FUZ_FACTOR = 2;
    return (p4 - (p4 % FUZ_FACTOR)) * 100000000 + (p3 - (p3 % FUZ_FACTOR))
    * 100000 + (p2 - (p2 % FUZ_FACTOR)) * 100
    + (p1 - (p1 % FUZ_FACTOR));
}
- (void)startProcessFile3 {
    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"Brad-Sucks--Total-Breakdown"
                                                              ofType:@"mp3"];
    NSError *outError;
    AVAsset *someAsset = [[AVURLAsset alloc] initWithURL:[NSURL fileURLWithPath:audioFilePath] options:nil];
    AVAssetReader *assetReader = [AVAssetReader assetReaderWithAsset:someAsset error:&outError];
    BOOL success = (assetReader != nil);
    AVAsset *localAsset = assetReader.asset;
    // Get the audio track to read.
    AVAssetTrack *audioTrack = [[localAsset tracksWithMediaType:AVMediaTypeAudio] objectAtIndex:0];
    // Decompression settings for Linear PCM
    NSDictionary *decompressionAudioSettings = @{ AVFormatIDKey : [NSNumber numberWithUnsignedInt:kAudioFormatLinearPCM] };
    // Create the output with the audio track and decompression settings.
    AVAssetReaderOutput *trackOutput = [AVAssetReaderTrackOutput assetReaderTrackOutputWithTrack:audioTrack outputSettings:decompressionAudioSettings];
    if ([assetReader canAddOutput:trackOutput])
        [assetReader addOutput:trackOutput];
    [assetReader startReading];
    BOOL done = NO;
    NSMutableArray *ray = [NSMutableArray new];
    // 8192
    int numberOfTotal = 0;
    while (!done)
    {
        numberOfTotal ++;
        // Copy the next sample buffer from the reader output.
        CMSampleBufferRef sampleBuffer = [trackOutput copyNextSampleBuffer];
        if (sampleBuffer)
        {
            
            CMItemCount numSamples = CMSampleBufferGetNumSamples(sampleBuffer);
//            printf(" %ld", numSamples);// 8192
            CMBlockBufferRef audioBuffer = CMSampleBufferGetDataBuffer(sampleBuffer);
            size_t lengthAtOffset;
            size_t totalLength;
            char *inSamples;
            CMBlockBufferGetDataPointer(audioBuffer, 0, &lengthAtOffset, &totalLength, &inSamples);
            CMAudioFormatDescriptionRef format = CMSampleBufferGetFormatDescription(sampleBuffer);
            const AudioStreamBasicDescription *desc = CMAudioFormatDescriptionGetStreamBasicDescription(format);
            assert(desc->mFormatID == kAudioFormatLinearPCM);
            
            if (desc->mChannelsPerFrame == 2 && desc->mBitsPerChannel == 16)
            {
                if (*inSamples == '\0')
                {
                    continue;
                }
                
                // Convert samples to floats
                float *samples = malloc(numSamples * sizeof(float));
                vDSP_vflt16((short *)inSamples, 1, samples, 1, numSamples);
                
                // Setup the FFT
                int fftRadix = log2(numSamples);
                int halfSamples = (int)(numSamples / 2);
                FFTSetup setup = vDSP_create_fftsetup(fftRadix, FFT_RADIX2);
                
                // Populate *window with the values for a hamming window function
                float *window = (float *)malloc(sizeof(float) * numSamples);
                vDSP_hamm_window(window, numSamples, 0);
                
                // Window the samples
                vDSP_vmul(samples, 1, window, 1, samples, 1, numSamples);
                
                // Define complex buffer
                COMPLEX_SPLIT A;
                A.realp = (float *) malloc(halfSamples * sizeof(float));
                A.imagp = (float *) malloc(halfSamples * sizeof(float));
                
                // Pack samples:
                vDSP_ctoz((COMPLEX*)samples, 2, &A, 1, numSamples/2);
                
                // Perform a forward FFT using fftSetup and A
                // Results are returned in A
                vDSP_fft_zrip(setup, &A, 1, fftRadix, FFT_FORWARD);
                
                // calculating square of magnitude for each value
                vDSP_zvmags(&A, 1, A.realp, 1, numSamples/2);
                
                // Inverse FFT
//                vDSP_fft_zrip(setup, &A, 1, fftRadix, kFFTDirection_Inverse);
//                
//                vDSP_ztoc(&A, 1, (COMPLEX *)samples, 2, numSamples/2);
//                
//                vDSP_Length lastZeroCrosssing;
//                vDSP_Length zeroCrossingCount;
//                vDSP_nzcros(samples, 1, numSamples, &lastZeroCrosssing, &zeroCrossingCount, numSamples);
                
                // Convert COMPLEX_SPLIT A result to magnitudes
                float amp[numSamples];
                amp[0] = A.realp[0]/(numSamples*2);
                
                // Find the max
                int maxIndex = 0;
                float maxMag = 0.0;
                
                // We can't detect anyting reliably above the Nyquist frequency
                // which is bin n / 2 and bin 0 should always empty.
                for(int i=1; i<halfSamples; i++)
                {
                    amp[i]=A.realp[i]*A.realp[i]+A.imagp[i]*A.imagp[i];
//                    printf("%f ",amp[i]);
//                    all = [NSString stringWithFormat:@"%@ %f",all, amp[i]];
                    [ray addObject:[NSNumber numberWithFloat: amp[i]]];
                    if (amp[i] > maxMag)
                    {
                        maxMag = amp[i];
                        maxIndex = i;
                    }
                }
//                NSLog(@"leng %lu", sizeof(amp)); //32768 if full
//                NSLog(@"max %f", maxMag); //32768 if full
                
            }
            
            
            CFRelease(sampleBuffer);
            sampleBuffer = NULL;
        }
        else
        {
            // Find out why the asset reader output couldn't copy another sample buffer.
            if (assetReader.status == AVAssetReaderStatusFailed)
            {
                NSError *failureError = assetReader.error;
                // Handle the error here.
            }
            else
            {
                // The asset reader output has read all of its samples.
                done = YES;
            }
        }
    }
    
  
    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
    NSString *documentsDirectory = [paths objectAtIndex:0];
    NSString *filePath = [NSString stringWithFormat:@"%@/plot.csv",documentsDirectory];
    [ray writeToFile:filePath atomically:YES];
    NSLog(@"numberOfTotal %d", numberOfTotal);
    NSLog(@"%lu", (unsigned long)ray.count);
}

- (void)ff:(float *)array {
    
}

- (void)startProcessFile2:(CAPAudioPlayer *) audioPlayer {
    
    ExtAudioFileRef audioFile;
    OSStatus status = noErr;
    // Open file
    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"Brad-Sucks--Total-Breakdown"
                                                              ofType:@"mp3"];
    NSURL *audioFileURL = [NSURL fileURLWithPath:audioFilePath];
    ExtAudioFileOpenURL((__bridge CFURLRef)audioFileURL, &audioFile);
    
    // Get files information
    AudioStreamBasicDescription fileAudioDescription;
    UInt32 size = sizeof(fileAudioDescription);
    ExtAudioFileGetProperty(audioFile,
                            kExtAudioFileProperty_FileDataFormat,
                            &size,
                            &fileAudioDescription);
    
    // Apply audio format
    ExtAudioFileSetProperty(audioFile,
                            kExtAudioFileProperty_ClientDataFormat,
                            sizeof(CAPAudioDescription),
                            &CAPAudioDescription);
    
    // Determine length in frames (in original file's sample rate)
    SInt64 fileLengthInFrames;
    size = sizeof(fileLengthInFrames);
    ExtAudioFileGetProperty(audioFile,
                            kExtAudioFileProperty_FileLengthFrames,
                            &size,
                            &fileLengthInFrames);
    NSLog(@"%lld", fileLengthInFrames);
    // Calculate the true length in frames, given the original and target sample rates
    fileLengthInFrames = ceil(fileLengthInFrames * (CAPAudioDescription.mSampleRate / fileAudioDescription.mSampleRate));
    
    // Prepare AudioBufferList: Interleaved data uses just one buffer, non-interleaved has two
    int numberOfBuffers = CAPAudioDescription.mFormatFlags & kAudioFormatFlagIsNonInterleaved ? CAPAudioDescription.mChannelsPerFrame : 1;
    int channelsPerBuffer = CAPAudioDescription.mFormatFlags & kAudioFormatFlagIsNonInterleaved ? 1 : CAPAudioDescription.mChannelsPerFrame;
    int bytesPerBuffer = CAPAudioDescription.mBytesPerFrame * (int)fileLengthInFrames;
    
    
    AudioBufferList *bufferList = malloc(sizeof(AudioBufferList) + (numberOfBuffers-1)*sizeof(AudioBuffer));
    
    bufferList->mNumberBuffers = numberOfBuffers;
    for ( int i=0; i<numberOfBuffers; i++ ) {
        if ( bytesPerBuffer > 0 ) {
            bufferList->mBuffers[i].mData = calloc(bytesPerBuffer, 1);
            if ( !bufferList->mBuffers[i].mData ) {
                for ( int j=0; j<i; j++ ) {
                    free(bufferList->mBuffers[j].mData);
                }
                free(bufferList);
                ExtAudioFileDispose(audioFile);
                printf("Could not allocate memory for buffer");
                exit(1);
                return;
            }
        } else {
            bufferList->mBuffers[i].mData = NULL;
        }
        bufferList->mBuffers[i].mDataByteSize = bytesPerBuffer;
        bufferList->mBuffers[i].mNumberChannels = channelsPerBuffer;
    }
    
    // Create a stack copy of the given audio buffer list and offset mData pointers, with offset in bytes
    char scratchBufferList_bytes[sizeof(AudioBufferList)+(sizeof(AudioBuffer)*(bufferList->mNumberBuffers-1))];
    memcpy(scratchBufferList_bytes, bufferList, sizeof(scratchBufferList_bytes));
    AudioBufferList * scratchBufferList = (AudioBufferList*)scratchBufferList_bytes;
    for ( int i=0; i<scratchBufferList->mNumberBuffers; i++ ) {
        scratchBufferList->mBuffers[i].mData = (char*)scratchBufferList->mBuffers[i].mData;
    }
    
    // Perform read in multiple small chunks (otherwise ExtAudioFileRead crashes when performing sample rate conversion)
    UInt32 readFrames = 0;
    float *audioTotal = malloc(fileLengthInFrames * sizeof(float));
//    NSMutableArray *ray = [NSMutableArray new];
    int index = 0;
    while (readFrames < fileLengthInFrames) {
        UInt32 framesLeftToRead = (UInt32)fileLengthInFrames - readFrames;
        UInt32 framesToRead = (16384 < framesLeftToRead) ? 16384 : framesLeftToRead;
        scratchBufferList->mNumberBuffers = bufferList->mNumberBuffers;
        scratchBufferList->mBuffers[0].mNumberChannels = bufferList->mBuffers[0].mNumberChannels;
        scratchBufferList->mBuffers[0].mData = bufferList->mBuffers[0].mData + (readFrames * CAPAudioDescription.mBytesPerFrame);
        scratchBufferList->mBuffers[0].mDataByteSize = framesToRead * CAPAudioDescription.mBytesPerFrame;

        status = ExtAudioFileRead(audioFile, &framesToRead, scratchBufferList);
        
        //extra
        SInt16 *inputFrames = (SInt16*)scratchBufferList->mBuffers[0].mData;
        SInt16 *convertedAudioData = (SInt16 *) malloc (sizeof(SInt16) * fileLengthInFrames);
        int x = 0;
        for (int i = 0; i < framesToRead * 2; i ++) {
            if(i%2 == 1){
                convertedAudioData[x] = inputFrames[i];
                x ++;
            }
        }
        
        int numSamples = (int)framesToRead;
        for(int i = 0; i < numSamples; ++i) {
            audioTotal[index] = (float)convertedAudioData[i];
            index ++;
        }
        free(convertedAudioData);
//
//        // Setup the length
//        vDSP_Length log2n = log2f(numSamples);
//        
//        // Calculate the weights array. This is a one-off operation.
//        FFTSetup fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
//        
//        // For an FFT, numSamples must be a power of 2, i.e. is always even
//        int nOver2 = numSamples/2;
//        
//        // Populate *window with the values for a hamming window function
//        float *window = (float *)malloc(sizeof(float) * numSamples);
//        vDSP_hamm_window(window, numSamples, 0);
//        // Window the samples
//        vDSP_vmul(samples, 1, window, 1, samples, 1, numSamples);
//        
//        // Define complex buffer
//        COMPLEX_SPLIT A;
//        A.realp = (float *) malloc(nOver2*sizeof(float));
//        A.imagp = (float *) malloc(nOver2*sizeof(float));
//        
//        // Pack samples:
//        // C(re) -> A[n], C(im) -> A[n+1]
//        vDSP_ctoz((COMPLEX*)samples, 2, &A, 1, numSamples/2);
//        
//        
//        //Perform a forward FFT using fftSetup and A
//        //Results are returned in A
//        vDSP_fft_zrip(fftSetup, &A, 1, log2n, FFT_FORWARD);
//        
//        //Convert COMPLEX_SPLIT A result to magnitudes
//        float amp[numSamples];
//        amp[0] = A.realp[0]/(numSamples*2);
//        
//        for(int i=1; i<numSamples; i++) {
//            amp[i]=sqrt(A.realp[i]*A.realp[i]+A.imagp[i]*A.imagp[i]);
//            [ray addObject:[NSNumber numberWithInt:amp[i]]];
////            printf("%f ",amp[i]);
//        }
        if ( framesToRead == 0 )
            break;
        
        readFrames += framesToRead;
    }
    
 
    free(bufferList);
//    NSArray *paths = NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES);
//    NSString *documentsDirectory = [paths objectAtIndex:0];
//    NSString *filePath = [NSString stringWithFormat:@"%@/plot.csv",documentsDirectory];
//    [ray writeToFile:filePath atomically:YES];
    // Clean up
    ExtAudioFileDispose(audioFile);

    // BufferList and readFrames are the audio we loaded
    audioPlayer->bufferList = bufferList;
    audioPlayer->frames = readFrames;
    

    //java
    int chunkSize = 4096;
    int totalSize = index;
    int sampledChunkSize = totalSize/chunkSize;
    float **result = malloc(sizeof(float) *sampledChunkSize);
    for (int i = 0; i < sampledChunkSize; i++) {
        result[i] = malloc(chunkSize * sizeof(float));
    }
    
    for(int j = 0;j < sampledChunkSize; j++) {
        float *complexArray = malloc(sizeof(float) * chunkSize);
        for(int i = 0; i < chunkSize; i++) {
            complexArray[i] = audioTotal[(j*chunkSize)+i];
        }

        // Setup the length
        vDSP_Length log2n = log2f(chunkSize);

        // Calculate the weights array. This is a one-off operation.
        FFTSetup fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);

        // For an FFT, numSamples must be a power of 2, i.e. is always even
        int nOver2 = chunkSize/2;

//        // Populate *window with the values for a hamming window function
//        float *window = (float *)malloc(sizeof(float) * chunkSize);
//        vDSP_hamm_window(window, chunkSize, 0);
//        // Window the samples
//        vDSP_vmul(complexArray, 1, window, 1, complexArray, 1, chunkSize);

        // Define complex buffer
        COMPLEX_SPLIT A;
        A.realp = (float *) malloc(nOver2*sizeof(float));
        A.imagp = (float *) malloc(nOver2*sizeof(float));

        // Pack samples:
        // C(re) -> A[n], C(im) -> A[n+1]
        vDSP_ctoz((COMPLEX*)complexArray, 2, &A, 1, chunkSize/2);
        
        //Perform a forward FFT using fftSetup and A
        //Results are returned in A
        vDSP_fft_zrip(fftSetup, &A, 1, log2n, FFT_FORWARD);
        
        //Convert COMPLEX_SPLIT A result to magnitudes
        float amp[chunkSize];
        amp[0] = A.realp[0]/(chunkSize*2);

        for(int i = 1; i < chunkSize; i++) {
            amp[i] = sqrt(A.realp[i]*A.realp[i] + A.imagp[i]*A.imagp[i]);
            printf("%f ",amp[i]);
            result[j][i] = amp[i];
        }
    }
    
    // result is complex matrix obtained in previous step
//    int points[sampledChunkSize][5];
//    for (int t = 0; t < sampledChunkSize; t++) {
//        for (int freq = 40; freq < 300 ; freq++) {
//            // Get the magnitude:
//            double mag = result[t][freq];
//            
//            // Find out which range we are in:
//            int indexFreq = [self getIndex:freq];
//            
//            // Save the highest magnitude and corresponding frequency:
//            if (mag > highscores[t][index]) {
//                points[t][indexFreq] = freq;
//            }
//        }
//        
//        // form hash tag
//        long h = hash(points[t][0], points[t][1], points[t][2], points[t][3]);
//    }

    

}

// find out in which range is frequency
- (int) getIndex:(int) freq {
    int RANGE [] = { 40, 80, 120, 180, 300 };
    int i = 0;
    while (RANGE[i] < freq)
        i++;
    return i;
}

- (void)startProcessFile {
    OSStatus err = noErr;
    AudioStreamBasicDescription fileFormat;
    NSString* audioFilePath = [[NSBundle mainBundle] pathForResource:@"Brad-Sucks--Total-Breakdown"
                                                     ofType:@"wav"];
    NSURL *audioFileURL = [NSURL fileURLWithPath:audioFilePath];
    AudioFileID afid;
//    OSStatus openAudioFileResult = AudioFileOpenURL((__bridge CFURLRef)audioFileURL, kAudioFileReadPermission, 0, &afid);
//
//    if (0 != openAudioFileResult) {
//        NSLog(@"An error occurred when attempting to open the audio file %@: %d", audioFilePath, (int)openAudioFileResult);
//    }
//
//    // get file size
//    UInt64 audioDataByteCount = 0;
//    UInt32 propertySize = sizeof(audioDataByteCount);
//    OSStatus getSizeResult = AudioFileGetProperty(afid, kAudioFilePropertyDataFormat, &propertySize, &audioDataByteCount);
//    
//    if (0 != getSizeResult)
//    {
//        NSLog(@"An error occurred when attempting to determine the size of audio file %@: %d", audioFilePath, (int)getSizeResult);
//    }
//
//    // get song sample rate
//    UInt32 propertySampleRateSize = sizeof(fileFormat);
//    err = AudioFileGetProperty(afid, kAudioFilePropertyDataFormat, &propertySampleRateSize, &fileFormat);
//    Float64 sampleRate = fileFormat.mSampleRate;
//    UInt32 channel = fileFormat.mChannelsPerFrame;
//    UInt32 bytesRead = (UInt32)audioDataByteCount;
//    void *audioData = malloc(bytesRead);
//    OSStatus readBytesResult = AudioFileReadBytes(afid, false, 0, &bytesRead, audioData);
//    
//    if (0 != readBytesResult) {
//        NSLog(@"An error occurred when attempting to read data from audio file %@: %d", audioFilePath, (int)readBytesResult);
//    }
//    if (audioData != NULL) {
//        NSLog(@"audio data pointer not empty");
//        AudioFileClose(afid);
//        free (audioData);
//    }
    
//    AudioFileClose(afid);
    NSData *fileData = [NSData dataWithContentsOfFile:audioFilePath];
    NSLog(@"fileData %@", [self sha1:fileData]);
    NSUInteger length = [fileData length];
    int *cdata = (int*)malloc(length);
    [fileData getBytes:(void*)cdata length:length];
    free(cdata);

    
    
    UInt32 thePropertySize = sizeof(fileFormat);
    ExtAudioFileRef extRef;
    err = ExtAudioFileOpenURL((__bridge CFURLRef _Nonnull)(audioFileURL), &extRef);
    
    // get channel and frame rate
    err = ExtAudioFileGetProperty(extRef, kExtAudioFileProperty_FileDataFormat, &thePropertySize, &fileFormat);
    BOOL interleaved = !(fileFormat.mFormatFlags & kAudioFormatFlagIsNonInterleaved);
    UInt32 channels = fileFormat.mChannelsPerFrame;
    Float64 sampleRate = fileFormat.mSampleRate;
    // get totaol frames
    UInt64 theFileLengthInFrames = 0;
    thePropertySize = sizeof(UInt64);
    err = ExtAudioFileGetProperty(extRef, kExtAudioFileProperty_FileLengthFrames, &thePropertySize, &theFileLengthInFrames);

    NSTimeInterval duration = theFileLengthInFrames/sampleRate;
    
//    UInt32 dataSize = theFileLengthInFrames) * theOutputFormat.mBytesPerFrame;

    float **data = (float **)malloc( sizeof(float*) * channels );
    for (int i = 0; i < channels; i++)
    {
        data[i] = (float *)malloc( sizeof(float) * 1024 );
    }
    
    SInt64 framesPerBuffer = ((SInt64) theFileLengthInFrames / 1024);
    SInt64 framesPerChannel = framesPerBuffer / channels;
    
    AudioBufferList *audioBufferList =  [self audioBufferListWithNumberOfFrames:(UInt32)framesPerBuffer
                                                                           numberOfChannels:channels
                                                                                interleaved:interleaved];
    
    
    // read through file and calculate rms at each point
    for (SInt64 i = 0; i < 1024; i++)
    {
        UInt32 bufferSize = (UInt32) framesPerBuffer;
        [self checkResult:ExtAudioFileRead(extRef,
                                                       &bufferSize,
                                                       audioBufferList)
                            operation:"Failed to read audio data from file waveform"];

        if (interleaved)
        {
            float *buffer = (float *)audioBufferList->mBuffers[0].mData;
            for (int channel = 0; channel < channels; channel++)
            {
                float channelData[framesPerChannel];
                for (int frame = 0; frame < framesPerChannel; frame++)
                {
                    channelData[frame] = buffer[frame * channels + channel];
                    NSLog(@"%f", channelData[frame] );
                }
                float rms = [self RMS:channelData length:(UInt32)framesPerChannel];
                data[channel][i] = rms;
            }
        }
        else
        {
            for (int channel = 0; channel < channels; channel++)
            {
                float *channelData = audioBufferList->mBuffers[channel].mData;
                float rms = [self RMS:channelData length:bufferSize];
                data[channel][i] = rms;
            }
        }
    }
}

- (void)checkResult:(OSStatus)result operation:(const char *)operation
{
    if (result == noErr) return;
    char errorString[20];
    // see if it appears to be a 4-char-code
    *(UInt32 *)(errorString + 1) = CFSwapInt32HostToBig(result);
    if (isprint(errorString[1]) && isprint(errorString[2]) && isprint(errorString[3]) && isprint(errorString[4]))
    {
        errorString[0] = errorString[5] = '\'';
        errorString[6] = '\0';
    } else
        // no, format it as an integer
        sprintf(errorString, "%d", (int)result);
    fprintf(stderr, "Error: %s (%s)\n", operation, errorString);

}

- (float)RMS:(float *)buffer   length:(int)bufferSize
{
    float sum = 0.0;
    for(int i = 0; i < bufferSize; i++)
        sum += buffer[i] * buffer[i];
    return sqrtf( sum / bufferSize);
}

- (AudioBufferList *)audioBufferListWithNumberOfFrames:(UInt32)frames
                                      numberOfChannels:(UInt32)channels
                                           interleaved:(BOOL)interleaved
{
    unsigned nBuffers;
    unsigned bufferSize;
    unsigned channelsPerBuffer;
    if (interleaved)
    {
        nBuffers = 1;
        bufferSize = sizeof(float) * frames * channels;
        channelsPerBuffer = channels;
    }
    else
    {
        nBuffers = channels;
        bufferSize = sizeof(float) * frames;
        channelsPerBuffer = 1;
    }
    
    AudioBufferList *audioBufferList = (AudioBufferList *)malloc(sizeof(AudioBufferList) + sizeof(AudioBuffer) * (channels-1));
    audioBufferList->mNumberBuffers = nBuffers;
    for(unsigned i = 0; i < nBuffers; i++)
    {
        audioBufferList->mBuffers[i].mNumberChannels = channelsPerBuffer;
        audioBufferList->mBuffers[i].mDataByteSize = bufferSize;
        audioBufferList->mBuffers[i].mData = calloc(bufferSize, 1);
    }
    return audioBufferList;
}

- (NSString *)sha1:(NSData *)data
{

    uint8_t digest[CC_SHA1_DIGEST_LENGTH];
    
    CC_SHA1(data.bytes, (CC_LONG)data.length, digest);
    
    NSMutableString *output = [NSMutableString stringWithCapacity:CC_SHA1_DIGEST_LENGTH * 2];
    
    for (int i = 0; i < CC_SHA1_DIGEST_LENGTH; i++)
    {
        [output appendFormat:@"%02x", digest[i]];
    }
    
    return output.uppercaseString;
}

#pragma mark - audio output functions -

void CAPStartAudioOutput (CAPAudioOutput *player) {
    OSStatus status = noErr;
    
    // Description for the output AudioComponent
    AudioComponentDescription outputcd = {
        .componentType = kAudioUnitType_Output,
        .componentSubType = kAudioUnitSubType_RemoteIO,
        .componentManufacturer = kAudioUnitManufacturer_Apple,
        .componentFlags = 0,
        .componentFlagsMask = 0
    };
    
    // Get the output AudioComponent
    AudioComponent comp = AudioComponentFindNext (NULL, &outputcd);
    if (comp == NULL) {
        printf ("can't get output unit");
        exit (-1);
    }
    
    // Create a new instance of the AudioComponent = the AudioUnit
    status = AudioComponentInstanceNew(comp, &player->outputUnit);
    CheckError (status, "Couldn't open component for outputUnit");
    
    
    // Set the stream format
    status = AudioUnitSetProperty(player->outputUnit,
                                  kAudioUnitProperty_StreamFormat,
                                  kAudioUnitScope_Input,
                                  0,
                                  &CAPAudioDescription,
                                  sizeof(CAPAudioDescription));
    CheckError (status,"kAudioUnitProperty_StreamFormat");
    
    
    // Set the render callback
    AURenderCallbackStruct input = {
        .inputProc = CAPRenderProc,
        .inputProcRefCon = player
    };
    
    status = AudioUnitSetProperty(player->outputUnit,
                                  kAudioUnitProperty_SetRenderCallback,
                                  kAudioUnitScope_Global,
                                  0,
                                  &input,
                                  sizeof(input));
    CheckError (status, "Could not set render callback");
    
    
    // Set the maximum frames per slice (not necessary)
    UInt32 framesPerSlice = 4096;
    status = AudioUnitSetProperty(player->outputUnit,
                                  kAudioUnitProperty_MaximumFramesPerSlice,
                                  kAudioUnitScope_Global,
                                  0,
                                  &framesPerSlice,
                                  sizeof(framesPerSlice));
    CheckError (status, "AudioUnitSetProperty(kAudioUnitProperty_MaximumFramesPerSlice");
    
    
    // Initialize the Audio Unit
    status = AudioUnitInitialize(player->outputUnit);
    CheckError (status, "Couldn't initialize output unit");
    
    
    // Start the Audio Unit (sound begins)
    status = AudioOutputUnitStart(player->outputUnit);
    CheckError (status, "Couldn't start output unit");
}

static OSStatus CAPRenderProc(void *inRefCon,
                              AudioUnitRenderActionFlags *ioActionFlags,
                              const AudioTimeStamp *inTimeStamp,
                              UInt32 inBusNumber,
                              UInt32 inNumberFrames,
                              AudioBufferList * ioData) {
    CAPAudioOutput *audioOutput = (CAPAudioOutput*)inRefCon;
    CAPAudioPlayer *audioPlayer = &audioOutput->player;
    
    UInt32 currentFrame = audioPlayer->currentFrame;
    UInt32 maxFrames = audioPlayer->frames;
    
    Float32 *outputData = (Float32*)ioData->mBuffers[0].mData;

    Float32 *inputData = (Float32*)audioPlayer->bufferList->mBuffers[0].mData;

        for (UInt32 frame = 0; frame < inNumberFrames; ++frame) {
            UInt32 outSample = frame * 2;
            UInt32 inSample = currentFrame * 2;
            
            (outputData)[outSample] = (inputData)[inSample];
            (outputData)[outSample+1] = (inputData)[inSample + 1];
            
            currentFrame++;
            currentFrame = currentFrame % maxFrames; // loop
    }
    
    audioPlayer->currentFrame = currentFrame;
    
    return noErr;
}

// generic error handler - if err is nonzero, prints error message and exits program.
static void CheckError(OSStatus error, const char *operation) {
    if (error == noErr) return;
    
    char str[20];
    // see if it appears to be a 4-char-code
    *(UInt32 *)(str + 1) = CFSwapInt32HostToBig(error);
    if (isprint(str[1]) && isprint(str[2]) && isprint(str[3]) && isprint(str[4])) {
        str[0] = str[5] = '\'';
        str[6] = '\0';
    } else
        // no, format it as an integer
        sprintf(str, "%d", (int)error);
    
    fprintf(stderr, "Error: %s (%s)\n", operation, str);
    
    exit(1);
}

void CAPDisposeAudioOutput(CAPAudioOutput *output) {
    AudioOutputUnitStop(output->outputUnit);
    AudioUnitUninitialize(output->outputUnit);
    AudioComponentInstanceDispose(output->outputUnit);
}

@end
