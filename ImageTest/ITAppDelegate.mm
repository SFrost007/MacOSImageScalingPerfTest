//
//  ITAppDelegate.m
//  ImageTest
//
//  Created by Simon Frost on 30/06/2014.
//  Copyright (c) 2014 Orangeninja. All rights reserved.
//

#import "ITAppDelegate.h"
#include <stddef.h>
#include <algorithm>
#include <cmath>

namespace util {
  typedef unsigned int uint32_t;
}


@implementation ITAppDelegate

typedef enum {
  Raw,
  NearestNeighbour,
  Bilinear,
  Fractional,
  Contributory,
  Enterprise,
  SCALING_ALGORITHM_MAX,
  Bicubic // Disabled
} ScalingAlgorithm;

struct Colour {
  Colour() : r(0), g(0), b(0), a(0) {}
  uint r,g,b,a;
};

struct entPoint {
  int x;
  int y;
};

struct ScalerCoords {
  entPoint srcSz, dstSz;
};

typedef uint32_t Pixel;
#define COLOUR_TYPE Colour
#define INPIXEL Pixel
#define OUTPIXEL Pixel
//#define HAS_ALPHA

#pragma mark - Application lifecycle

- (NSString*) nameForScalingAlgorithm:(ScalingAlgorithm)algorithm
{
  switch (algorithm) {
    case Raw:
      return @"Raw";
    case NearestNeighbour:
      return @"Nearest";
    case Bilinear:
      return @"Bilinear";
    case Bicubic:
      return @"Bicubic";
    case Contributory:
      return @"Contributory";
    case Fractional:
      return @"Fractional";
    case Enterprise:
      return @"Enterprise";
    case SCALING_ALGORITHM_MAX:
      assert(0);
  }
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
  self.window.delegate = self;
  
  // Update the algorithm switcher items
  [self.algorithmSwitch setSegmentCount:SCALING_ALGORITHM_MAX];
  CGFloat segmentWidth = CGRectGetWidth(self.algorithmSwitch.frame) / SCALING_ALGORITHM_MAX - 2.0;
  for (int i=0; i<SCALING_ALGORITHM_MAX; i++) {
    NSString *name = [self nameForScalingAlgorithm:(ScalingAlgorithm)i];
    [self.algorithmSwitch setWidth:segmentWidth forSegment:i];
    [self.algorithmSwitch setLabel:name forSegment:i];
  }
  self.algorithmSwitch.selectedSegment = Fractional;
  
  // Render the startup image
  [self rerenderImage];
}

- (BOOL)applicationShouldTerminateAfterLastWindowClosed:(NSApplication *)sender
{
  return YES;
}

#pragma mark - UI handling

- (void) windowDidResize:(NSNotification *)notification
{
  [self rerenderImage];
}

- (IBAction)algorithmChanged:(id)sender
{
  [self rerenderImage];
}

- (void)rerenderImage
{
  ScalingAlgorithm algorithm = (ScalingAlgorithm)[self.algorithmSwitch selectedSegment];
  
  NSDate *startTime = [NSDate date];
  NSBitmapImageRep *imageRep = [self imageRepForScalingAlgorithm:algorithm];
  NSTimeInterval duration = [[NSDate date] timeIntervalSinceDate:startTime];
  
  NSString *timeString = [NSString stringWithFormat:@"Compute time: %f", duration];
  [self.computeTimeLabel setStringValue:timeString];
  NSString *resString = [NSString stringWithFormat:@"Resolution: (%d,%d)",
                         (int)self.imageView.frame.size.width, (int)self.imageView.frame.size.height];
  [self.resolutionLabel setStringValue:resString];
  
  NSImage *image = [[NSImage alloc] init];
  [image addRepresentation:imageRep];
  self.imageView.image = image;
}

#pragma mark - Scaling algorithms

- (NSBitmapImageRep*) imageRepForScalingAlgorithm:(ScalingAlgorithm)algorithm
{
  NSString *path = [[NSBundle mainBundle] pathForResource:@"small" ofType:@"png"];
  NSData *inputData = [NSData dataWithContentsOfFile:path];
  NSBitmapImageRep *inputImage = [[NSBitmapImageRep alloc] initWithData:inputData];
  
  const int targetWidth = self.imageView.frame.size.width;
  const int targetHeight = self.imageView.frame.size.height;
  
  NSBitmapImageRep *outputImage = [[NSBitmapImageRep alloc] initWithBitmapDataPlanes:nil
                                                                          pixelsWide:targetWidth
                                                                          pixelsHigh:targetHeight
                                                                       bitsPerSample:inputImage.bitsPerSample
                                                                     samplesPerPixel:inputImage.samplesPerPixel
                                                                            hasAlpha:inputImage.hasAlpha
                                                                            isPlanar:inputImage.isPlanar
                                                                      colorSpaceName:inputImage.colorSpaceName
                                                                        bitmapFormat:inputImage.bitmapFormat
                                                                         bytesPerRow:(inputImage.bitsPerPixel * targetWidth) / 8
                                                                        bitsPerPixel:inputImage.bitsPerPixel];
  
  
  
  const std::size_t srcWidth = (int)inputImage.pixelsWide;
  const std::size_t srcHeight = (int)inputImage.pixelsHigh;
  const std::size_t destWidth = targetWidth;
  const std::size_t destHeight = targetHeight;
  
  const std::size_t srcStridePixels = srcWidth;
  const std::size_t destStridePixels = destWidth;
  
  const double scaleX = (double)destWidth / (double)srcWidth;
  const double scaleY = (double)destHeight / (double)srcHeight;
  const double scaleFactor = std::min(scaleX, scaleY);
  
  std::size_t numRows = std::floor(srcHeight * scaleFactor);
  std::size_t numCols = std::floor(srcWidth * scaleFactor);
  
  std::size_t marginX = 0, marginY = 0;
  if (scaleX < scaleY) {
    marginY = (destHeight - (scaleFactor * srcHeight)) / 2;
  } else {
    marginX = (destWidth - (scaleFactor * srcWidth)) / 2;
  }
  
  const std::size_t destLengthBytes = destStridePixels * sizeof(Pixel) * destHeight;
  const std::size_t destLengthPixels = destLengthBytes / sizeof(Pixel);
  
  const Pixel* srcBegin = reinterpret_cast<uint32_t*>(inputImage.bitmapData);
  Pixel* destBegin = reinterpret_cast<uint32_t*>(outputImage.bitmapData);
  Pixel* destEnd = destBegin + destLengthPixels;
  
  // Windows 8 requires us to fill the margins, otherwise window ends up transparent
  // TODO: Is this too expensive? We only really need to write in the margins
  Pixel* wipePixel = destBegin;
  while (wipePixel < destEnd) {
    *(wipePixel++) = mask_op(0);
  }
  destBegin += (marginY * destStridePixels);
  
  if (algorithm == Raw) {
    outputImage = inputImage;
    destBegin = destEnd;
    
  } else if (algorithm == NearestNeighbour) {
    for (std::size_t row = 0; row < numRows; ++row) {
      const std::size_t srcRow = std::floor(row / scaleFactor);
      const Pixel* srcRowBegin = srcBegin + (srcRow * srcStridePixels);
      Pixel* destPixel = destBegin + (row * destStridePixels) + marginX;
      
      for (std::size_t col=0; col<numCols; ++col, ++destPixel) {
        const std::size_t srcCol = std::floor(col / scaleFactor);
        const Pixel* srcPixel = srcRowBegin + srcCol;
        
        *destPixel = mask_op(*srcPixel);
      }
    }
    
  } else if (algorithm == Bilinear) {
    util::uint32_t p0, p1, p2, p3, x, y;
    std::size_t index;
    double x_ratio = ((double)(srcWidth-1)) / numCols;
    double y_ratio = ((double)(srcHeight-1)) / numRows;
    double x_diff, y_diff;
    int r, g, b;

    for (std::size_t row=0; row<numRows; row++) {
      Pixel* destPixel = destBegin + (row * destStridePixels) + marginX;
      for (std::size_t col=0; col<numCols; col++, destPixel++) {
        x = (int)(x_ratio * col);
        y = (int)(y_ratio * row);
        x_diff = (x_ratio * col) - x;
        y_diff = (y_ratio * row) - y;
        index = (y*srcWidth+x);
        p0 = *(srcBegin + index);
        p1 = *(srcBegin + index+1);
        p2 = *(srcBegin + index+srcWidth);
        p3 = *(srcBegin + index+srcWidth+1);
        
        b = PixelB(p0)*(1-x_diff)*(1-y_diff)
          + PixelB(p1)*(x_diff)*(1-y_diff)
          + PixelB(p2)*(y_diff)*(1-x_diff)
          + PixelB(p3)*(x_diff*y_diff);
        
        g = PixelG(p0)*(1-x_diff)*(1-y_diff)
          + PixelG(p1)*(x_diff)*(1-y_diff)
          + PixelG(p2)*(y_diff)*(1-x_diff)
          + PixelG(p3)*(x_diff*y_diff);
        
        r = PixelR(p0)*(1-x_diff)*(1-y_diff)
          + PixelR(p1)*(x_diff)*(1-y_diff)
          + PixelR(p2)*(y_diff)*(1-x_diff)
          + PixelR(p3)*(x_diff*y_diff);
        
        Pixel out = ((r<<16)&0xff0000) | ((g<<8)&0xff00) | (b&0xff);
        *destPixel = mask_op(out);
      }
    }
    
  } else if (algorithm == Bicubic) {
    
    /*
    
    std::size_t a, b, c, index;
    Pixel C[5];
    Pixel Cc, d0, d2, d3, a0, a1, a2, a3;
    int x, y;
    double dx, dy;
    
    double x_ratio = (double) srcWidth / numCols;
    double y_ratio = (double) destWidth / numRows;
    
    for (std::size_t row=0; row<numRows; row++) {
      for (std::size_t col=0; col<numCols; col++) {
        x = (int)(x_ratio*col);
        y = (int)(y_ratio*row);
        
        dx = x_ratio * col - x;
        dy = y_ratio * row - y;
        
        index = y * srcWidth + x;
        a = y * srcWidth + (x+1);
        b = (y+1) * srcWidth + x;
        c = (y+1) * srcWidth + (x+1);
        
        for (std::size_t k=0; k<3; k++) {
          for (std::size_t jj=0; jj<=3; jj++) {
            d0 = srcBegin[(y-1+jj)*srcWidth + (x-1) + k] - srcBegin[(y-1+jj)*srcWidth + x + k];
            d2 = srcBegin[(y-1+jj)*srcWidth + (x+1) + k] - srcBegin[(y-1+jj)*srcWidth + x + k];
            d3 = srcBegin[(y-1+jj)*srcWidth + (x+2) + k] - srcBegin[(y-1+jj)*srcWidth + x + k];
            a0 = srcBegin[(y-1+jj)*srcWidth + x + k]; // TODO: This can be optimised
            a1 = -1.0/3*d0 + d2 - 1.0/6*d3;
            a2 = 1.0/2*d0 - 1.0/2*d2;
            a3 = -1.0/6*d0 - 1.0/2*d2 + 1.0/6*d3;
            C[jj] = a0 + a1*dx + a2*dx*dx + a3*dx*dx*dx;
            
            d0 = C[0]-C[1];
            d2 = C[2]-C[1];
            d3 = C[3]-C[1];
            a0 = C[1];
            
            a1 = -1.0/3*d0 + d2 - 1.0/6*d3;
            a2 = 1.0/2*d0 + 1.0/2*d2;
            a3 = -1.0/6*d0 - 1.0/2*d2 + 1.0/6*d3;
            Cc = a0 + a1*dy + a2*dy*dy + a3*dy*dy*dy;
            destBegin[row*numCols + col + k] = Cc;
          }
        }
      }
    }
    
    /*
    
    util::uint32_t a, b, c, d, x, y;
    int index;
    double x_ratio = ((double)(srcWidth-1)) / numCols;
    double y_ratio = ((double)(srcHeight-1)) / numRows;
    double x_diff, y_diff;
    
    for (std::size_t row=1; row<numRows-3; row++) {
      Pixel *destPixel = destBegin + (row * destStridePixels) + marginX;
      for (std::size_t col=1; col<numCols-3; col++, destPixel++) {
        x = x_ratio * col;
        y = y_ratio * row;
        x_diff = x_ratio * col - x;
        y_diff = y_ratio * row - y;
        
        index = fmax((y-1),0)*(int)srcWidth;
        const Pixel* p = (srcBegin + index + x);
        a = interpolatePixel(*(p-1), *p, *(p+1), *(p+2), x_diff);
        
        index += srcWidth;
        p = (srcBegin + index + x);
        b = interpolatePixel(*(p-1), *p, *(p+1), *(p+2), x_diff);
        
        index += srcWidth;
        p = (srcBegin + index + x);
        c = interpolatePixel(*(p-1), *p, *(p+1), *(p+2), x_diff);
        
        index += srcWidth;
        p = (srcBegin + index + x);
        d = interpolatePixel(*(p-1), *p, *(p+1), *(p+2), x_diff);
        
        *destPixel = mask_op(interpolatePixel(a, b, c, d, y_diff));
      }
    }
    
    //*/
  
  } else if (algorithm == Contributory) {
    const double pixelRatio = 1 / scaleFactor;
    double realY, realX;
    int startX, startY, endX, endY;
    int r, g, b;
    for (std::size_t row=0; row<numRows; ++row) {
      Pixel* destPixel = destBegin + (row * destStridePixels) + marginX;
      
      realY = (double)row * pixelRatio;
      startY = std::floor(realY);
      endY = std::min(std::floor(realY + pixelRatio), (double)srcWidth-1);
      for (std::size_t col=0; col<=numCols; ++col, ++destPixel) {
        
        realX = (double)col * pixelRatio;
        startX = std::floor(realX);
        endX = std::min(std::floor(realX + pixelRatio), (double)srcWidth-1);
        const double numPixels = (endX+1-startX)*(endY+1-startY);
        
        r = g = b = 0;
        
        for (int srcY=startY; srcY<=endY; ++srcY) {
          const Pixel* srcRowBegin = srcBegin + (srcY * srcStridePixels);
          for (int srcX=startX; srcX<=endX; srcX++) {
            const Pixel srcPixel = *(srcRowBegin + srcX);
            b += (srcPixel&0xff);
            g += ((srcPixel>>8)&0xff);
            r += ((srcPixel>>16)&0xff);
          }
        }
        b /= numPixels;
        g /= numPixels;
        r /= numPixels;
        
        Pixel out = ((r<<16)&0xff0000) | ((g<<8)&0xff00) | (b&0xff);
        *destPixel = mask_op(out);
      }
    }
    
    
  } else if (algorithm == Fractional) {
    const double pixelRatio = 1 / scaleFactor;
    double realY, realX, xFrac, yFrac;
    int startX, startY, endX, endY;
    uint32_t r, g, b;
    for (std::size_t row=0; row<numRows; ++row) {
      Pixel* destPixel = destBegin + (row * destStridePixels) + marginX;
      
      realY = (double)row * pixelRatio;
      startY = std::floor(realY);
      endY = std::min(std::floor(realY + pixelRatio), (double)srcHeight-1);
      for (std::size_t col=0; col<=numCols; ++col, ++destPixel) {
        r = g = b = 0;
        double totalOverlap = 0;
        
        realX = (double)col * pixelRatio;
        startX = std::floor(realX);
        endX = std::min(std::floor(realX + pixelRatio), (double)srcWidth-1);
        for (int srcY=startY; srcY<=endY; ++srcY) {
          if (srcY == startY) {
            yFrac = fmod(srcY+1, pixelRatio);
          } else if (srcY == endY) {
            yFrac = 1-fmod(srcY+1, pixelRatio);
          } else {
            yFrac = 1;
          }
          yFrac = std::max(0.0, yFrac);
          
          const Pixel* srcRowBegin = srcBegin + (srcY * srcStridePixels);
          for (int srcX=startX; srcX<=endX; srcX++) {
            if (srcX == startX) {
              xFrac = fmod(srcX+1, pixelRatio);
            } else if (srcX == endX) {
              xFrac = 1-fmod(srcX+1, pixelRatio);
            } else {
              xFrac = 1;
            }
            xFrac = std::max(0.0, xFrac);
            
            const double overlap = xFrac * yFrac;
            const Pixel srcPixel = *(srcRowBegin + srcX);
            b += (srcPixel&0xff)*overlap;
            g += ((srcPixel>>8)&0xff)*overlap;
            r += ((srcPixel>>16)&0xff)*overlap;
            totalOverlap += overlap;
          }
        }
        b /= totalOverlap;
        g /= totalOverlap;
        r /= totalOverlap;
        
        Pixel outP = ((r<<16)&0xff0000) | ((g<<8)&0xff00) | (b&0xff);
        *destPixel = mask_op(outP);
      }
    }
  
  } else if (algorithm == Enterprise) {
    
    // Set up some variables to match names in enterprise algorithm
    Colour* rowBuf_ = new Colour[numCols];
    Pixel *dstPels_ = destBegin + marginX;
    const Pixel *srcPels_ = srcBegin;
    const std::size_t dstStride = destStridePixels;
    const std::size_t srcStride = srcStridePixels;
    ScalerCoords coords_;
    coords_.srcSz.x = (int)srcWidth;
    coords_.srcSz.y = (int)srcHeight;
    coords_.dstSz.x = (int)numCols;
    coords_.dstSz.y = (int)numRows;
	
#ifdef HAS_ALPHA
    COLOUR_TYPE* rowBuf = (COLOUR_TYPE*)rowBuf_;
#else
    COLOUR_TYPE* rowBuf = rowBuf_;
#endif
    
    OUTPIXEL* dstPels = (OUTPIXEL*)dstPels_;
    const INPIXEL* srcPels = (INPIXEL*)srcPels_;
    int srcNormY = 0;
    int dstNormY = 0;
    COLOUR_TYPE* rowEnd = rowBuf + numCols;
    
    // Pre-calculate multipliers to avoid having to use division per-pixel
    int colourMulX = (1<<24) / coords_.srcSz.x;
    int colourMulY = (1<<24) / coords_.srcSz.y;
    int extraX = coords_.srcSz.x / 2;
    int extraY = coords_.srcSz.y / 2;
    
    int nextSrcNormY = srcNormY + coords_.dstSz.y;
    int nextDstNormY = dstNormY + coords_.srcSz.y;
    
    for (std::size_t y=0; y<numRows; y++) {
      // Loop accumulating source-row contents
      while (true) {
        const INPIXEL* srcPos = srcPels;
        int srcNormX = 0;
        int dstNormX = 0;
        
        for (COLOUR_TYPE* rowPos=rowBuf; rowPos<rowEnd; rowPos++) {
          COLOUR_TYPE dstRgb;
          int nextDstNormX = dstNormX+coords_.srcSz.x;
          
          // Loop accumulating source-row contents
          while (true) {
            int nextSrcNormX = srcNormX + coords_.dstSz.x;
            COLOUR_TYPE srcRgb;
            PEL2RGB(*srcPos, &srcRgb);
            int overlap_ = overlap(srcNormX, nextSrcNormX,
                                   dstNormX, nextDstNormX);
            dstRgb.r += srcRgb.r * overlap_;
            dstRgb.g += srcRgb.g * overlap_;
            dstRgb.b += srcRgb.b * overlap_;
#ifdef HAS_ALPHA
            dstRgb.a += srcRgb.a * overlap_;
#endif
            // Exit loop when current source pixel's right edge is beyond
            // the right edge of the destination pixel
            if (nextSrcNormX >= nextDstNormX)
              break;
            srcNormX = nextSrcNormX;
            srcPos++;
          }
          int overlap_ = overlap(srcNormY, nextSrcNormY,
                                 dstNormY, nextDstNormY);
          rowPos->r += (((dstRgb.r + extraX) * colourMulX) >> 24) * overlap_;
          rowPos->g += (((dstRgb.g + extraX) * colourMulX) >> 24) * overlap_;
          rowPos->b += (((dstRgb.b + extraX) * colourMulX) >> 24) * overlap_;
#ifdef HAS_ALPHA
          rowPos->a += (((dstRgb.a + extraX) * colourMulX) >> 24) * overlap_;
#endif
          dstNormX = nextDstNormX;
        }
        
        // Exit loop when bottom of source row is further down
        // the screen than the bottom of the destination row
        if (nextSrcNormY >= nextDstNormY)
          break;
        
        // Move to the next row of source pixels
        srcNormY = nextSrcNormY;
        nextSrcNormY += coords_.dstSz.y;
        srcPels += srcStride;
      }
      
      // Apply the accumulated colour values to the destination pixels and
      // clear the row-accumulation buffer
      OUTPIXEL *dstPos = dstPels;
      COLOUR_TYPE c;
      for (COLOUR_TYPE* rowPos=rowBuf; rowPos<rowEnd; rowPos++) {
        c.r = ((rowPos->r + extraY) * colourMulY) >> 24;
        c.g = ((rowPos->g + extraY) * colourMulY) >> 24;
        c.b = ((rowPos->b + extraY) * colourMulY) >> 24;
#ifdef HAS_ALPHA
        c.a = ((rowPos->a + extraY) * colourMulY) >> 24;
#endif
        
        *(dstPos++) = RGB2PEL(c);
        rowPos->r = rowPos->g = rowPos->b = 0;
#ifdef HAS_ALPHA
        rowPos->a = 0;
#endif
      }
      
      // Move on to the next destination row
      dstPels += dstStride;
      dstNormY = nextDstNormY;
      nextDstNormY += coords_.srcSz.y;
    }
    
    delete [] rowBuf;
  } // End algorithm if
  
  return outputImage;
}


#pragma mark - Utility functions

Pixel mask_op(const Pixel pixel) {
  return pixel | 0; // TODO: Incorrect
}

inline static int overlap(int srcXY, int srcXYright, int dstXY, int dstXYright)
{
  int tl = (dstXY >= srcXY) ? dstXY : srcXY;
  int br = (dstXYright >= srcXYright) ? srcXYright : dstXYright;
  return br - tl;
}

uint interpolateValue(uint v0, uint v1, uint v2, uint v3, float t)
{
  int p = (v3 - v2) - (v0 - v1);
  int q = (v0 - v1) - p;
  int r = v2 - v0;
  int s = v1;
  float tSqrd = t * t;
  return (p * (tSqrd * t)) + (q * tSqrd) + (r * t) + s;
}

Pixel interpolatePixel(Pixel p0, Pixel p1, Pixel p2, Pixel p3, float t)
{
  return interpolateValue(PixelR(p0), PixelR(p1), PixelR(p2), PixelR(p3), t) << 16 |
         interpolateValue(PixelG(p0), PixelG(p1), PixelG(p2), PixelG(p3), t) << 8 |
         interpolateValue(PixelB(p0), PixelB(p1), PixelB(p2), PixelB(p3), t);
}

uint PixelR(Pixel p)
{
  return (p>>16)&0xff;
}

uint PixelG(Pixel p)
{
  return (p>>8)&0xff;
}

uint PixelB(Pixel p)
{
  return p&0xff;
}

void PEL2RGB(const Pixel inPixel, Colour *outRgb)
{
  outRgb->r = PixelR(inPixel);
  outRgb->g = PixelG(inPixel);
  outRgb->b = PixelB(inPixel);
}

Pixel RGB2PEL(const COLOUR_TYPE c)
{
  return ((c.r<<16)&0xff0000) | ((c.g<<8)&0xff00) | (c.b&0xff);
}

@end
