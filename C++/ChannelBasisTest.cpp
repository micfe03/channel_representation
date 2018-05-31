//
//  ChannelBasisTest.cpp
//  ChannelBasisTest
//
//  Created by Michael Felsberg on 9/14/11.
//  Copyright 2012 Link�ping University. All rights reserved.
//

#include "ChannelBasis.hpp"
#include <iostream>
#include <fstream>
#include <time.h>

int main (int argc, const char * argv[])
{

//	int axdim = 7;
//	int Nch = 16;
//	cvl::CombinedChannelBasis jointEncoder;
//    //cvl::BsplineChannelBasis  jointEncoder1D;
//    cvl::Cos2ChannelBasis     jointEncoder1D;
//    jointEncoder1D.setParameters(Nch,-M_PI,M_PI,1); // bounded, not modular!
//    std::vector<cvl::ChannelBasis*> jointEncChs(axdim, &jointEncoder1D); // 7
//    jointEncoder.setParameters(jointEncChs);
//    jointEncoder.displayParameters();
//
//    std::vector<float> K(axdim);
//    K[0] = 0.3f;
//    K[1] = 0.31f;
//    K[2] = 0.32f;
//    K[3] = 0.33f;
//    K[4] = 0.34f;
//    K[5] = 0.35f;
//    K[6] = 0.36f;
//
//    cv::Mat_<float> chCoefficients(jointEncoder.getNrChannels(),1);
//
//    jointEncoder.encode(chCoefficients, K);
////    std::ifstream matdump("/Users/micha/Desktop/outmatbsp_10.dump", std::ios::binary);
//    std::streamsize s = jointEncoder.getNrChannels();		// 16^7
////    matdump.read((char *)chCoefficients.data, s*sizeof(float));
//
//    std::ofstream matdump("/Users/micha/Desktop/outmatbsp_10_mfe.dump", std::ios::binary);
//    const char* sdata = (char*) chCoefficients.data;
//    matdump.write(sdata, s*sizeof(float));
//
//    cv::Mat_<float> res;
//
//    jointEncoder.decode(res, chCoefficients.t());
//
//    std::cout << "Decoded " << res << '\n';
//
//    return 0;


    // random generation of parameters
    srand((unsigned int)time(0));
    int nMax = 30;
    int nCCB = (rand() % (nMax-2)) + 3;
    int nBCB = (rand() % (nMax-2)) + 3;
    float rCMax = (float)rand()/(float)RAND_MAX;
    float rCMin = -(float)rand()/(float)RAND_MAX;
    float rBMax = (float)rand()/(float)RAND_MAX;
    float rBMin = -(float)rand()/(float)RAND_MAX;

    cv::Mat_<float> x(2,2), y(2,2), r(2,2);
    x(0,0) = rCMin+(rCMax-rCMin)*(float)rand()/(float)RAND_MAX;
    y(0,0) = rBMin+(rBMax-rBMin)*(float)rand()/(float)RAND_MAX;
    r(0,0) = 1.0f;
    x(0,1) = rCMin+(rCMax-rCMin)*(float)rand()/(float)RAND_MAX;
    y(0,1) = rBMin+(rBMax-rBMin)*(float)rand()/(float)RAND_MAX;
    r(0,1) = 0.9f;
    x(1,0) = rCMin+(rCMax-rCMin)*(float)rand()/(float)RAND_MAX;
    y(1,0) = rBMin+(rBMax-rBMin)*(float)rand()/(float)RAND_MAX;
    r(1,0) = 0.8f;
    x(1,1) = rCMin+(rCMax-rCMin)*(float)rand()/(float)RAND_MAX;
    y(1,1) = rBMin+(rBMax-rBMin)*(float)rand()/(float)RAND_MAX;
    r(1,1) = 0.7f;

    cv::Mat_<cv::Vec<float, 2> > xr(2,2), yr(2,2), xy(2,2);
    cv::Mat_<cv::Vec<float, 3> > xyr(2,2);

    cv::Mat XR[] = {x, r};
    cv::Mat YR[] = {y, r};
    cv::Mat XY[] = {x, y};
    cv::Mat XYR[] = {x, y, r};

    cv::merge(XR, 2, xr);
    cv::merge(YR, 2, yr);
    cv::merge(XY, 2, xy);
    cv::merge(XYR, 3, xyr);

    cv::Mat_<float> result;

    std::cout << "ChannelBasisTest\n\n";
    
    std::cout << "Cos2-channels with " << nCCB << " coefficients from " << rCMin << " to " << rCMax << "\n";

    cvl::Cos2ChannelBasis CCB;
    CCB.setParameters(nCCB, rCMin, rCMax);
    CCB.displayParameters();
    std::cout << "Norm is " << CCB.getNorm() << "\n\n";
    
    std::cout << "Bspline-channels with " << nBCB << " coefficients from " << rBMin << " to " << rBMax << "\n";    
    
    cvl::BsplineChannelBasis BCB;
    BCB.setParameters(nBCB, rBMin, rBMax);
    BCB.displayParameters();
    std::cout << "Norm is " << BCB.getNorm() << "\n\n";
    
    std::cout << "Cross product of channels\n";
    
    cvl::CombinedChannelBasis CBCB;
    std::vector<cvl::ChannelBasis*> chBV(2);
    chBV[0] = &CCB;
    chBV[1] = &BCB;
    CBCB.setParameters(chBV);
    CBCB.displayParameters();
    std::cout << "Norm is " << CBCB.getNorm() << "\n\n";

    std::cout << "Channel vectors (non-sparse)\n\n";
    
    std::cout << "Encode and decode " << x(cv::Range(0,1),cv::Range(0,1)) << " as cos2 channel\n";
    cvl::ChannelVector cCV(&CCB);
    cCV.addSample(x(cv::Range(0,1),cv::Range(0,1)));
    std::cout << "Vector is " << cCV << "\n";
    std::flush(std::cout);
    cCV.decode(result);
    std::cout << "Decoding gives "<< result << "\n\n";
    result.setTo(0.f);

    std::cout << "Add " << xr(cv::Range(0,1),cv::Range(1,2)) << " as cos2 channel and decode\n";
    cCV.addSample(xr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << cCV << "\n";
    cCV.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);
        
    cCV.setTo(0.0f);
    cCV.setChannelBasis(&CCB);
    
    std::cout << "Encode and decode " << xr(cv::Range(0,1),cv::Range(0,2)) << " as cos2 channel\n";
    cCV.addSample(xr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << cCV << "\n";
    cCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);
    
    cvl::ChannelVector cCVsum(&CCB,(cCV(cv::Range(0,1),cv::Range::all())+cCV(cv::Range(1,2),cv::Range::all())));
    std::cout << "Average the two channels gives " << cCVsum << "\n";
    cCVsum.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Encode and decode " << y(cv::Range(0,1),cv::Range(0,1)) << " as bspline channel\n";
    cvl::ChannelVector bCV(&BCB);
    bCV.addSample(y(cv::Range(0,1),cv::Range(0,1)));
    std::cout << "Vector is " << bCV << "\n";
    std::flush(std::cout);
    bCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);
    
    std::cout << "Add " << yr(cv::Range(0,1),cv::Range(1,2)) << " as bspline channel and decode\n";
    bCV.addSample(yr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << bCV << "\n";
    bCV.decode(result, 2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    bCV.setTo(0.0f);

    std::cout << "Encode and decode " << yr(cv::Range(0,1),cv::Range(0,2)) << " as bspline channel\n";
    bCV.addSample(yr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << bCV << "\n";
    bCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    cvl::ChannelVector bCVsum(&BCB,(bCV(cv::Range(0,1),cv::Range::all())+bCV(cv::Range(1,2),cv::Range::all())));
    std::cout << "Average the two channels gives " << bCVsum << "\n";
    bCVsum.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Encode and decode " << xy(cv::Range(0,1),cv::Range(0,1)) << " as combined channels\n";
    cvl::ChannelVector cbCV(&CBCB);
    cbCV.addSample(xy(cv::Range(0,1),cv::Range(0,1)));
    std::cout << "Vector is " << cbCV << "\n";
    cbCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Add " << xyr(cv::Range(0,1),cv::Range(1,2)) << " as combined channel and decode\n";
    cbCV.addSample(xyr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << cbCV << "\n";
    cbCV.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    cbCV.setTo(0.0f);

    std::cout << "Encode and decode " << xyr(cv::Range(0,1),cv::Range(0,2)) << " as combined channel\n";
    cbCV.addSample(xyr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << cbCV << "\n";
    cbCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    cvl::ChannelVector cbCVsum(&CBCB,(cbCV(cv::Range(0,1),cv::Range::all())+cbCV(cv::Range(1,2),cv::Range::all())));
    std::cout << "Average the two channels gives " << cbCVsum << "\n";
    cbCVsum.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);
    
    std::cout << "Channel vectors (sparse)\n\n";
    
    std::cout << "Encode and decode " << x(cv::Range(0,1),cv::Range(0,1)) << " as cos2 channel\n";
    cvl::ChannelSVector cSCV(&CCB);
    cSCV.addSample(x(cv::Range(0,1),cv::Range(0,1)));
    std::cout << "Vector is " << cSCV << "\n";
    cSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Add " << xr(cv::Range(0,1),cv::Range(1,2)) << " as cos2 channel and decode\n";
    cSCV.addSample(xr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << cSCV << "\n";
    cSCV.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    cSCV.clear();
    
    std::cout << "Encode and decode " << xr(cv::Range(0,1),cv::Range(0,2)) << " as cos2 channel\n";
    cSCV.addSample(xr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << cSCV << "\n";
    cSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

   // sum skipped as sparse things cannot be accessed by column
    
    std::cout << "Encode and decode " << y(cv::Range(0,1),cv::Range(0,1)) << " as bspline channel\n";
    cvl::ChannelSVector bSCV(&BCB);
    bSCV.addSample(y(cv::Range(0,1),cv::Range(0,1)));
    std::cout << "Vector is " << bSCV << "\n";
    bSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Add " << yr(cv::Range(0,1),cv::Range(1,2)) << " as bspline channel and decode\n";
    bSCV.addSample(yr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << bSCV << "\n";
    bSCV.decode(result, 2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);
    
    bSCV.clear();

    std::cout << "Encode and decode " << yr(cv::Range(0,1),cv::Range(0,2)) << " as bspline channel\n";
    bSCV.addSample(yr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << bSCV << "\n";
    bSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Encode and decode " << xy(cv::Range(0,1),cv::Range(0,1))  << " as combined channels\n";
    cvl::ChannelSVector cbSCV(&CBCB);
    cbSCV.addSample(xy(cv::Range(0,1),cv::Range(0,1)) );
    std::cout << "Vector is " << cbSCV << "\n";
    cbSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    std::cout << "Add " << xyr(cv::Range(0,1),cv::Range(1,2)) << " as combined channel and decode\n";
    cbSCV.addSample(xyr(cv::Range(0,1),cv::Range(1,2)));
    std::cout << "Vector is " << cbSCV << "\n";
    cbSCV.decode(result,2);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    cbSCV.clear();

    std::cout << "Encode and decode " << xyr(cv::Range(0,1),cv::Range(0,2)) << " as combined channel\n";
    cbSCV.addSample(xyr(cv::Range(0,1),cv::Range(0,2)));
    std::cout << "Vector is " << cbSCV << "\n";
    cbSCV.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";

    cv::Mat image;
    image = cv::imread("house_orig.png");
    CV_Assert(~image.empty());
    image.pop_back(10); // just to check not to mix x and y
    cv::cvtColor(image,image,cv::COLOR_RGB2GRAY);
    cv::namedWindow("Original Image");
    cv::imshow("Original Image", image);
//    std::cout << image.rows << " , " << image.cols << " , " << image.channels() << " , " << (int)image.at<uchar>(100,100) << std::endl;
    cv::Mat_<float> fImage(image.size());
    cv::MatConstIterator_<uchar> it=image.begin<uchar>(), it_end=image.end<uchar>();
    cv::MatIterator_<float> it2=fImage.begin();
    for ( ; it!=it_end; ++it, ++it2)
    	*it2=((float)*it)/255.0f;
    cv::imshow("Original Image", fImage);
    std::cout << fImage.rows << " , " << fImage.cols << " , " << fImage.channels() << " , " << fImage(100,250) << std::endl;
    cvl::Cos2ChannelBasis CCBB;
    CCBB.setParameters(10, 0, 1);
    CCBB.displayParameters();
    cvl::ChannelVector cCVB(&CCBB);
//    waitKey(0);
    cCVB.addSample(fImage);
    cv::namedWindow("Channel Image");
    //for filtering, but GaussianBlur must do some internal split!!
    //cv::Mat_<cv::Vec<float, 10> > channelImage=cCVB.reshape(10,fImage.rows);
    cv::Mat channelImage;
    cCVB.channelImage(channelImage);
    cv::GaussianBlur(channelImage, channelImage, cv::Size(7,7), 1.5, 1.5);
    //cv::Mat_<float> tmpCoeffs = channelImage.reshape(1,fImage.rows*fImage.cols);
    //cvl::ChannelVector cCVBs(&CCBB,tmpCoeffs);
    cv::Mat fImage2;
    cCVB.decode(fImage2,2);
    std::vector<cv::Mat_<float> > fImageLayers(2);
    cv::split(fImage2.reshape(0,fImage.rows),fImageLayers); //TODO: workaround because class cannot know the true dimensions (=bad!)
    cv::imshow("Channel Image", fImageLayers[0]);
    cv::waitKey(0);

    // CCFM test

    cvl::Cos2ChannelBasis CCBC; // column basis
    CCBC.setParameters(14, 0, fImage.cols-1);
    CCBC.displayParameters();

    cvl::Cos2ChannelBasis CCBR; // row basis
    CCBR.setParameters(12, 0, fImage.rows-1);
    CCBR.displayParameters();

    CCBB.setParameters(10, 0, 1);
    CCBB.displayParameters();

    cvl::CombinedChannelBasis CCFM;
    std::vector<cvl::ChannelBasis*> chBCCFM(3);
    chBCCFM[0] = &CCBB;
    chBCCFM[1] = &CCBR;
    chBCCFM[2] = &CCBC;
    CCFM.setParameters(chBCCFM);
    CCFM.displayParameters();
    std::cout << "Norm is " << CCFM.getNorm() << "\n\n";

    cvl::ChannelVector cCCFM(&CCFM);
    cCCFM.addSample(fImage);

//    std::cout << "Vector is " << cCCFM << "\n\n";

    cCCFM.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    // sparse CCFM

    fImage = cv::Mat::ones(11,13,CV_32F);
    fImage(5,6) = 0.f;

    CCBC.setParameters(3, 0, fImage.cols-1);
    CCBC.displayParameters();

    CCBR.setParameters(4, 0, fImage.rows-1);
    CCBR.displayParameters();

    CCBB.setParameters(5, 0, 1);
    CCBB.displayParameters();

    chBCCFM[0] = &CCBB;
    chBCCFM[1] = &CCBR;
    chBCCFM[2] = &CCBC;
    CCFM.setParameters(chBCCFM);
    CCFM.displayParameters();
    std::cout << "Norm is " << CCFM.getNorm() << "\n\n";

    cvl::ChannelSVector cSCCFM(&CCFM);
    cSCCFM.addSample(fImage);
    for (int xy = 0 ; xy < 3*4 ; xy++) {
        	cSCCFM.ref(3+5*xy,0) = 0.f;
        	cSCCFM.ref(4+5*xy,0) = 0.f;
    }
    std::cout << "Vector is " << cSCCFM << "\n";
    cSCCFM.decode(result);
    std::cout << "Decoding gives " << result << "\n\n";
    result.setTo(0.f);

    return 0;
}
