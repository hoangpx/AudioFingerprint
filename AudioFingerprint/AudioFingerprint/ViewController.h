//
//  ViewController.h
//  AudioFingerprint
//
//  Created by Hoang Pham on 1/11/17.
//  Copyright © 2017 Hoang Pham. All rights reserved.
//
#include <stdio.h>
#include <AudioToolbox/AudioToolbox.h>
#import <UIKit/UIKit.h>
typedef struct CAPAudioPlayer {
    AudioBufferList *bufferList;
    UInt32 frames;
    UInt32 currentFrame;
} CAPAudioPlayer;
typedef struct CAPAudioOutput
{
    AudioUnit outputUnit;
    double startingFrameCount;
    CAPAudioPlayer player;
} CAPAudioOutput;
@interface ViewController : UIViewController


@end

