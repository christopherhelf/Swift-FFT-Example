//
//  fft.swift
//  ffttest
//
//  Created by Christopher Helf on 17.08.15.
//  Copyright (c) 2015 Christopher Helf. All rights reserved.
//  Adapted From https://gerrybeauregard.wordpress.com/2013/01/28/using-apples-vdspaccelerate-fft/

import Foundation
import Accelerate

class FFT {
    
    private func getFrequencies(N: Int, fps: Double) -> [Double] {
    
        // Create an Array with the Frequencies
        let freqs = (0..<N/2).map {
            fps/Double(N)*Double($0)
        }
        
        return freqs
    }
    
    private func generateBandPassFilter(freqs: [Double]) -> ([Double], Int, Int) {
        
        var minIdx = freqs.count+1
        var maxIdx = -1
        
        let bandPassFilter: [Double] = map(enumerate(freqs)) { (index, element) in
            if (element >= self.lowerFreq && element <= self.higherFreq) {
                return 1.0
            } else {
                return 0.0
            }
        }
        
        for(var i=0; i<bandPassFilter.count; i++) {
            if(bandPassFilter[i]==1.0) {
                if(i<minIdx || minIdx == freqs.count+1) {
                    minIdx=i
                }
                if(i>maxIdx || maxIdx == -1) {
                    maxIdx=i
                }
            }
        }
        
        assert(maxIdx != -1)
        assert(minIdx != freqs.count+1)
        
        return (bandPassFilter, minIdx, maxIdx)
        
    }
    
    func calculate(_values: [Double], fps: Double) {
        
        // ----------------------------------------------------------------
        // Copy of our input
        // ----------------------------------------------------------------
        var values = _values
        
        // ----------------------------------------------------------------
        // Size Variables
        // ----------------------------------------------------------------
        let N = values.count
        let N2 = vDSP_Length(N/2)
        let LOG_N = vDSP_Length(log2(Float(values.count)))
        
        // ----------------------------------------------------------------
        // FFT & Variables Setup
        // ----------------------------------------------------------------
        let fftSetup: FFTSetupD = vDSP_create_fftsetupD(LOG_N, FFTRadix(kFFTRadix2))
        
        // We need complex buffers in two different formats!
        var tempComplex : [DSPDoubleComplex] = [DSPDoubleComplex](count: N/2, repeatedValue: DSPDoubleComplex())
        
        var tempSplitComplexReal : [Double] = [Double](count: N/2, repeatedValue: 0.0)
        var tempSplitComplexImag : [Double] = [Double](count: N/2, repeatedValue: 0.0)
        var tempSplitComplex : DSPDoubleSplitComplex = DSPDoubleSplitComplex(realp: &tempSplitComplexReal, imagp: &tempSplitComplexImag)
        
        // For polar coordinates
        var mag : [Double] = [Double](count: N/2, repeatedValue: 0.0)
        var phase : [Double] = [Double](count: N/2, repeatedValue: 0.0)
        
        // ----------------------------------------------------------------
        // Forward FFT
        // ----------------------------------------------------------------
        
        var valuesAsComplex : UnsafeMutablePointer<DSPDoubleComplex>?
        values.withUnsafeBufferPointer { (resultPointer: UnsafeBufferPointer<Double>) -> Void in
            valuesAsComplex = UnsafeMutablePointer<DSPDoubleComplex>( resultPointer.baseAddress )
        }
        
        // Scramble-pack the real data into complex buffer in just the way that's
        // required by the real-to-complex FFT function that follows.
        vDSP_ctozD(valuesAsComplex!, 2, &tempSplitComplex, 1, N2);
        
        // Do real->complex forward FFT
        vDSP_fft_zripD(fftSetup, &tempSplitComplex, 1, LOG_N, FFTDirection(FFT_FORWARD));
        
        // ----------------------------------------------------------------
        // Get the Frequency Spectrum
        // ----------------------------------------------------------------
        
        var fftMagnitudes = [Double](count:N/2, repeatedValue:0.0)
        vDSP_zvmagsD(&tempSplitComplex, 1, &fftMagnitudes, 1, N2);
        
        // vDSP_zvmagsD returns squares of the FFT magnitudes, so take the root here
        let roots = sqrt(fftMagnitudes)
        
        // Normalize the Amplitudes
        var fullSpectrum = [Double](count:N/2, repeatedValue:0.0)
        vDSP_vsmulD(roots, vDSP_Stride(1), [1.0 / Double(N)], &fullSpectrum, 1, N2)
        
        // ----------------------------------------------------------------
        // Convert from complex/rectangular (real, imaginary) coordinates
        // to polar (magnitude and phase) coordinates.
        // ----------------------------------------------------------------
        
        vDSP_zvabsD(&tempSplitComplex, 1, &mag, 1, N2);
        
        // Beware: Outputted phase here between -PI and +PI
        // https://developer.apple.com/library/prerelease/ios/documentation/Accelerate/Reference/vDSPRef/index.html#//apple_ref/c/func/vDSP_zvphasD
        vDSP_zvphasD(&tempSplitComplex, 1, &phase, 1, N2);
        
        // Save this variable for output
        var fullPhases = phase
        
        // ----------------------------------------------------------------
        // Bandpass Filtering
        // ----------------------------------------------------------------
        
        // Get the Frequencies for the current Framerate
        let freqs = getFrequencies(N,fps: fps)
        // Get a Bandpass Filter
        let bandPassFilter = generateBandPassFilter(freqs)
        
        // Multiply phase and magnitude with the bandpass filter
        mag = mul(mag, y: bandPassFilter.0)
        phase = mul(phase, y: bandPassFilter.0)
        
        // Output Variables
        var filteredSpectrum = mul(fullSpectrum, y: bandPassFilter.0)
        var filteredPhase = phase
        
        // ----------------------------------------------------------------
        // Determine Maximum Frequency
        // ----------------------------------------------------------------
        var maxFrequencyResult = max(filteredSpectrum)
        var maxFrequency = freqs[maxFrequencyResult.1]
        var maxPhase = filteredPhase[maxFrequencyResult.1]
        
        println("Amplitude: \(maxFrequencyResult.0)")
        println("Frequency: \(maxFrequency)")
        println("Phase: \(maxPhase + M_PI/2)")
        
        // ----------------------------------------------------------------
        // Convert from polar coordinates back to rectangular coordinates.
        // ----------------------------------------------------------------
        
        tempSplitComplex = DSPDoubleSplitComplex(realp: &mag, imagp: &phase)
        
        var complexAsValue : UnsafeMutablePointer<Double>?
        tempComplex.withUnsafeBufferPointer { (resultPointer: UnsafeBufferPointer<DSPDoubleComplex>) -> Void in
            complexAsValue = UnsafeMutablePointer<Double>( resultPointer.baseAddress )
        }
        
        vDSP_ztocD(&tempSplitComplex, 1, &tempComplex, 2, N2);
        vDSP_rectD(complexAsValue!, 2, complexAsValue!, 2, N2);
        vDSP_ctozD(&tempComplex, 2, &tempSplitComplex, 1, N2);
        
        // ----------------------------------------------------------------
        // Do Inverse FFT
        // ----------------------------------------------------------------
        
        // Create result
        var result : [Double] = [Double](count: N, repeatedValue: 0.0)
        var resultAsComplex : UnsafeMutablePointer<DSPDoubleComplex>?

        result.withUnsafeBufferPointer { (resultPointer: UnsafeBufferPointer<Double>) -> Void in
            resultAsComplex = UnsafeMutablePointer<DSPDoubleComplex>( resultPointer.baseAddress )
        }

        // Do complex->real inverse FFT.
        vDSP_fft_zripD(fftSetup, &tempSplitComplex, 1, LOG_N, FFTDirection(FFT_INVERSE));
        
        // This leaves result in packed format. Here we unpack it into a real vector.
        vDSP_ztocD(&tempSplitComplex, 1, resultAsComplex!, 2, N2);
        
        // Neither the forward nor inverse FFT does any scaling. Here we compensate for that.
        var scale : Double = 0.5/Double(N);
        vDSP_vsmulD(&result, 1, &scale, &result, 1, vDSP_Length(N));
 
        // Print Result
        for(var k=0; k<N; k++) {
            println("\(k)   \(values[k])     \(result[k])")
        }

    }
    
    // The bandpass frequencies
    let lowerFreq : Double = 3
    let higherFreq: Double = 5
    
    // Some Math functions on Arrays
    func mul(x: [Double], y: [Double]) -> [Double] {
        var results = [Double](count: x.count, repeatedValue: 0.0)
        vDSP_vmulD(x, 1, y, 1, &results, 1, vDSP_Length(x.count))
        return results
    }
    
    func sqrt(x: [Double]) -> [Double] {
        var results = [Double](count:x.count, repeatedValue:0.0)
        vvsqrt(&results, x, [Int32(x.count)])
        return results
    }
    
    func max(x: [Double]) -> (Double, Int) {
        var result: Double = 0.0
        var idx : vDSP_Length = vDSP_Length(0)
        vDSP_maxviD(x, 1, &result, &idx, vDSP_Length(x.count))
        return (result, Int(idx))
    }
    
}
