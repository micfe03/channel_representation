#include "ChannelBasis.hpp"

//Constructors
cvl::ChannelBasis::ChannelBasis(){
	//Setting default values
	m_nChans  = 11;
	m_minV    = 0;
	m_maxV    = 0;
	m_modular = 0; 
	m_bounded = 1; 
	m_width   = (float)M_PI/3;
	//Update the mapping constants
	updateInputScaling();    
}

const int cvl::ChannelBasis::getNrChannels() const {
	return m_nChans;
}

const float cvl::ChannelBasis::getNorm() const {
	return m_norm;
}

const float cvl::ChannelBasis::getMaxV() const {
	return m_maxV;
}

const float cvl::ChannelBasis::getMinV() const {
	return m_minV;
}

void cvl::ChannelBasis::updateInputScaling(){
	m_scaling = (m_nChans-2)/(m_maxV-m_minV);
	m_offset  = 0.5;	
	return;
}

void cvl::ChannelBasis::scaleInput(float &val) const {
	if (m_bounded){
		val=std::max(val,m_minV);
		val=std::min(val,m_maxV);
	}
	val = m_scaling*(val-m_minV)+m_offset;
	return;
}

void cvl::ChannelBasis::setParameters(const int nChannels, const float minValue, const float maxValue, const int modularChannels, const int boundedChannels){
	// Setting parameters 
	m_nChans  = nChannels;
	m_minV    = minValue;
	m_maxV    = maxValue;
	m_modular = modularChannels; 
	m_bounded = boundedChannels; //bounded support or let the edge channels have infinite bandwidth
	updateInputScaling();    
	return;
}

void cvl::ChannelBasis::displayParameters() const {
	displayBasisType();
	std::cout << " nchans       : " << m_nChans << std::endl;
	std::cout << " [minv,maxv]  : [" << m_minV << "," << m_maxV << "]" << std::endl;
	std::cout << " modular      : " << m_modular << std::endl;
	std::cout << " bounded      : " << m_bounded << std::endl;
	return;
}

// --------------  Cos2ChannelBasis -------------------

cvl::Cos2ChannelBasis::Cos2ChannelBasis(){
	m_norm=1.5;
}

void cvl::Cos2ChannelBasis::encode(cv::Mat_<float> &chVec, const float cVal) const {
	//Map the input value into the domain spanned by the channels  
	float val = cVal;
	scaleInput(val);
	//Setting all channel values to zero
	chVec=cv::Mat_<float>::zeros(m_nChans,1);

	//The channel centers are at placed at integers [0,numberOfChannels-1]  	
	//So we find the index of the center activated channel and the distance to it	
	int i_centerCoeff = (int)(val+0.5f);	// round(val)
	i_centerCoeff=std::min(i_centerCoeff,m_nChans-2);
	float d = val -(float)i_centerCoeff;

	//We are in the middle of the interval spanned by the channels
	chVec(i_centerCoeff-1)   = std::pow(std::cos((1+d)*m_width),2.0f);
	chVec(i_centerCoeff)     = std::pow(std::cos(d*m_width),2.0f);
	chVec(i_centerCoeff+1)   = std::pow(std::cos((1-d)*m_width),2.0f);

	return;
}

void cvl::Cos2ChannelBasis::encode(cv::SparseMat_<float> &chSVec, const float cVal) const {
	//Map the input value into the domain spanned by the channels
	float val = cVal;
	scaleInput(val);
	int sizes[1] = {m_nChans};
	chSVec.create(1,sizes);

	//The channel centers are at placed at integers [0,numberOfChannels-1]  	
	//So we find the index of the center activated channel and the distance to it	
	int i_centerCoeff = (int)(val+0.5f);	// round(val)
	i_centerCoeff=std::min(i_centerCoeff,m_nChans-2);
	float d = val -(float)i_centerCoeff;

	//We are in the middle of the interval spanned by the channels
	// bugfix for OpenCVs bug on 1D sparse mats
	chSVec.ref(i_centerCoeff)     = std::pow(std::cos(d*m_width),2.0f);
	i_centerCoeff--;
	chSVec.ref(i_centerCoeff)   = std::pow(std::cos((1+d)*m_width),2.0f);
	i_centerCoeff+=2;
	chSVec.ref(i_centerCoeff)   = std::pow(std::cos((1-d)*m_width),2.0f);

	return;
}

void cvl::Cos2ChannelBasis::decode(cv::Mat_<float> &res, const cv::Mat_<float> &chCoeff, const int nModes) const {
	//Decoding of a cos2 channel representation
	CV_Assert(m_nChans>2);
	std::vector<float> resVec(m_nChans-2,0);
	std::vector<float> confVec(m_nChans-2,0);
	std::complex<float> imag(0, 1);
	std::vector<float> v_tmpConf(nModes,-1);

	//Temporary result and certainty to be of the right size
	std::vector<float> result(nModes,0);
	std::vector<float> certainty(nModes,-1);
	//Setting the certainty to -1

	//Setting res to correct size
	res.create(1,nModes*2);

	//Decode
	float argument = 0;
	for(int i=0;i<m_nChans-2;i++){
		argument = std::arg((std::complex<float>)chCoeff.at<float>(i)+((std::complex<float>)chCoeff.at<float>(i+1)*std::exp(imag*2.f*m_width))+((std::complex<float>)chCoeff.at<float>(i+2)*std::exp(imag*4.f*m_width)));
		if(argument<0)
			argument+=2.f*3.14159265358979323846f;

		resVec[i]  = (float)i+1.f/(2.f*m_width)*argument;
		confVec[i] = 1.f/3.f *(chCoeff.at<float>(i)+chCoeff.at<float>(i+1)+chCoeff.at<float>(i+2))*3.f;
	}

	//Checking that the result is in the valid range
	for(int i=0;i<m_nChans-2;i++){
		if((resVec[i]<i+0.5) || (resVec[i]>i+1.5)){
			resVec[i]  = 0;
			confVec[i] = 0;
		}
		else
			resVec[i]=(resVec[i]-m_offset)/m_scaling+m_minV;
		//Mapping back to [minV maxV]
	}

	sortModes(resVec, confVec, result, certainty); // todo: opencv has inbuilt sort method - use that?
	for(int i=0; i<nModes; i++){ // this could be argued to be more or less optimal,
		// but if  different number of modes are to be aligned, this is better!
		res(0,2*i) = result[i];
		res(0,2*i+1) = certainty[i]/m_norm;
	}
	return;
}

void cvl::sortModes(std::vector<float> &resVec, std::vector<float> &confVec, std::vector<float> &res, std::vector<float> &certainty)
{
	int nModes = (int)res.size();
	float tmpC;

	//Two different sorting loops depending on if we just need the first mode or more
	//If we need just one mode, no need to sort. We just identify to strongest mode
	if(nModes==1) 
	{ 
		for(unsigned i=0;i<resVec.size();i++){
			if(confVec[i]>certainty[0])
			{
				res[0]=resVec[i];
				certainty[0]=confVec[i];
			}
		}
		return;
	}

	//We need to do some sorting to find the n strongest modes
	for(unsigned i=0;i<resVec.size();i++)
	{
		if(confVec[i]>certainty[nModes-1])
		{
			certainty[nModes-1]=confVec[i];
			res[nModes-1]=resVec[i];
			{
				for(int j=nModes-1;j>0;j--)
				{
					if(certainty[j]<certainty[j-1])
						break;
					tmpC=certainty[j-1];
					certainty[j-1]=certainty[j];
					certainty[j]=tmpC;

					tmpC=res[j-1];
					res[j-1]=res[j];
					res[j]=tmpC;
				}
			}
		}
	}
	return;
}

void cvl::Cos2ChannelBasis::displayBasisType() const
{
	std::cout << " Channel type : cos2" << std::endl;  
	return;
};

// --------------  BsplineChannelBasis -------------------

cvl::BsplineChannelBasis::BsplineChannelBasis(){
	m_norm=1.0;
}

void cvl::BsplineChannelBasis::encode(cv::Mat_<float> &chVec, const float cVal) const {
	//Map the input value into the domain spanned by the channels
	float val = cVal;
	scaleInput(val);
	//Setting all channel values to zero
	chVec=cv::Mat_<float>::zeros(m_nChans,1);

	//The channel centers are at placed at integers [0,numberOfChannels-1]  	
	//So we find the index of the center activated channel and the distance to it	
	int i_centerCoeff = (int)(val+0.5f);	
	i_centerCoeff=std::min(i_centerCoeff,m_nChans-2);
	float d = val -(float)i_centerCoeff;

	//We are in the middle of the interval spanned by the channels
	chVec(i_centerCoeff-1)   = 0.5f*std::pow(0.5f-d,2);
	chVec(i_centerCoeff)     = 0.75f-std::pow(d,2);
	chVec(i_centerCoeff+1)   = 0.5f*std::pow(0.5f+d,2);

	return;
}

void cvl::BsplineChannelBasis::encode(cv::SparseMat_<float> &chSVec, const float cVal) const {
	//Map the input value into the domain spanned by the channels  
	float val = cVal;
	scaleInput(val);
	int sizes[1] = {m_nChans};
	chSVec.create(1,sizes);

	//The channel centers are at placed at integers [0,numberOfChannels-1]  	
	//So we find the index of the center activated channel and the distance to it	
	int i_centerCoeff = (int)(val+0.5f);	
	i_centerCoeff=std::min(i_centerCoeff,m_nChans-2);
	float d = val -(float)i_centerCoeff;

	//We are in the middle of the interval spanned by the channels
	// bugfix for OpenCVs bug on 1D sparse mats
	chSVec.ref(&i_centerCoeff)     = 0.75f-std::pow(d,2);
	i_centerCoeff--;
	chSVec.ref(&i_centerCoeff)   = 0.5f*std::pow(0.5f-d,2);
	i_centerCoeff+=2;
	chSVec.ref(&i_centerCoeff)   = 0.5f*std::pow(0.5f+d,2);

	return;
}

void cvl::BsplineChannelBasis::decode(cv::Mat_<float> &res, const cv::Mat_<float> &chCoeff, const int nModes) const {
	//Decoding of a Bspline channel representation
	CV_Assert(m_nChans>2);
	std::vector<float> resVec(m_nChans-2,0);
	std::vector<float> confVec(m_nChans-2,0);

	//Temporary result and certainty to be of the right size
	std::vector<float> result(nModes,0);
	std::vector<float> certainty(nModes,-1);
	//Setting the certainty to -1

	//Seeting res to correct size
	res.create(1,nModes*2);

	float n, p, q, m, S;
	const float z = 2.0f*std::sqrt(2.0f)-3.0f;
	std::vector<float> ch(m_nChans+2);

	// recursive filtering
	ch[0] = 0.0;
	for (int cc=1; cc<m_nChans+1; cc++)
		ch[cc] = chCoeff.at<float>(cc-1) + z*ch[cc-1];
	ch[m_nChans+1] = z*ch[m_nChans];
	ch[m_nChans+1] = 8.0f*z/(z*z-1.0f)*ch[m_nChans+1];
	for (int cc=m_nChans; cc>=0; cc--)
		ch[cc] = -8.0f*z*ch[cc] + z*ch[cc+1];
	// zero detection + choosing the one with highest energy
	for (int cc=2; cc<m_nChans; cc++) {
		n = ch[cc+2] - ch[cc-2] + 2.0f * ( ch[cc-1] - ch[cc+1] );
		p = ch[cc] - 0.5f * ( ch[cc-2] + ch[cc+2]);
		q = 0.25f * ( ch[cc+2] - ch[cc-2] ) + 1.5f * ( ch[cc+1] - ch[cc-1] );
		S = p * p - q * n;
		S = (float)( p - std::sqrt( (S>0.f)?S:0.f ) ) / ( (n!=0)?n:1.0e-15f );
		m = ( ch[cc-2] + ch[cc+2] ) / 48.0f + ( ch[cc-1] + ch[cc+1] ) / 2.0f + 23.0f / 24.0f * ch[cc];
		confVec[cc-2] = 1-24.0f/23.0f*((fabs(S)<0.505f)?( 23.0f/24.0f -  m - S * ( q + S * ( S * n / 3.0f - p) ) / 2.0f ):(23.0f/24.0f));
		resVec[cc-2] = ( S + cc - m_offset - 1.0f ) / m_scaling + m_minV;
	}

	sortModes(resVec, confVec, result, certainty); // todo: opencv has inbuilt sort method - use that?
	for(int i=0; i<nModes; i++){ // this could be argued to be more or less optimal,
		// but if  different number of modes are to be aligned, this is better!
		res(0,2*i) = result[i];
		res(0,2*i+1) = certainty[i];
	}
	return;
}

void cvl::BsplineChannelBasis::displayBasisType() const
{
	std::cout << " Channel type : bspline" << std::endl;  
	return;
};

// ---------------	Combined Channel Basis -------------

cvl::CombinedChannelBasis::CombinedChannelBasis(){}

void cvl::CombinedChannelBasis::setParameters(const std::vector<ChannelBasis*> chBasisVector) {
	m_chBasisVector = chBasisVector;
	m_nChans = 1;
	m_norm = 1;
	m_nChansVec.resize(m_chBasisVector.size());
	for (unsigned int k=0; k<m_chBasisVector.size(); k++) {
		const int nPartChannels = m_chBasisVector[k]->getNrChannels();
		CV_Assert(nPartChannels > 0);
		CV_Assert(m_nChans <= (std::numeric_limits<int>::max()) / nPartChannels); // if ( b > INT_MAX / a ) // a * b would overflow
		m_nChans *= nPartChannels;
		m_norm *= m_chBasisVector[k]->getNorm();
		m_nChansVec[k] = nPartChannels;
	}
	return;
}

void cvl::CombinedChannelBasis::getNrChannelsVec(std::vector<int> &nrChansVec) const {
	nrChansVec = m_nChansVec;
}

void cvl::CombinedChannelBasis::encode(cv::Mat_<float> &chCoefficients, const std::vector<float> &vals) const {
	if (vals.size()==m_chBasisVector.size()) {
		cvl::ChannelSVector tmpxChannel(m_chBasisVector[0]);
		cv::Mat_<float> mVals(1,1);
		mVals(0,0) = vals[0];
		tmpxChannel.addSample(mVals);
		cvl::ChannelSVector tmpOP(m_chBasisVector[0]);
		for (unsigned int k=1; k<vals.size(); k++) {
			cvl::ChannelSVector tmpyChannel(m_chBasisVector[k]);
			mVals(0,0) = vals[k];
			tmpyChannel.addSample(mVals);
			int nrXChannels = tmpxChannel.size(0);
			int sizes[1] = {nrXChannels*tmpyChannel.size(0)};
			tmpOP.create(1,sizes);
			for (cv::SparseMatIterator_<float> spr=tmpyChannel.begin(); spr!=tmpyChannel.end(); ++spr)
				for (cv::SparseMatIterator_<float> spc=tmpxChannel.begin(); spc!=tmpxChannel.end(); ++spc) {
					int idx = spr.node()->idx[0]*nrXChannels+spc.node()->idx[0];
					tmpOP.ref(&idx) = (*spr)*(*spc);
				}
			tmpxChannel = tmpOP;
		}
		tmpxChannel.normalize();
		tmpxChannel.copyTo(chCoefficients);
	}
	return;
}

void cvl::CombinedChannelBasis::encode(cv::SparseMat_<float> &chCoefficients, const std::vector<float> &vals) const {
	if (vals.size()==m_chBasisVector.size()) {
		cvl::ChannelSVector tmpxChannel(m_chBasisVector[0]);
		cv::Mat_<float> mVals(1,1);
		mVals(0,0) = vals[0];
		tmpxChannel.addSample(mVals);
		cvl::ChannelSVector tmpOP(m_chBasisVector[0]);
		for (unsigned int k=1; k<vals.size(); k++) {
			cvl::ChannelSVector tmpyChannel(m_chBasisVector[k]);
			mVals(0,0) = vals[k];
			tmpyChannel.addSample(mVals);
			int nrXChannels = tmpxChannel.size(0);
			int sizes[] = {nrXChannels*tmpyChannel.size(0)};
			tmpOP.create(1,sizes);
			for (cv::SparseMatIterator_<float> spr=tmpyChannel.begin(); spr!=tmpyChannel.end(); ++spr)
				for (cv::SparseMatIterator_<float> spc=tmpxChannel.begin(); spc!=tmpxChannel.end(); ++spc) {
					int dxs[1];
					dxs[0] = spr.node()->idx[0]*nrXChannels+spc.node()->idx[0];
					tmpOP.ref(dxs) = (*spr)*(*spc);
				}
			tmpxChannel = tmpOP;
		}
		tmpxChannel.normalize();
		tmpxChannel.copyTo(chCoefficients);
	}
	return;
}

void cvl::CombinedChannelBasis::displayBasisType() const {
	std::cout << "Combined Channel Basis" << std::endl;
}

void cvl::CombinedChannelBasis::decode(cv::Mat_<float> &res, const cv::Mat_<float> &chCoeff, const int nrModes) const {
	int nrXChannels;
	cv::Mat_<float> xChannel;
	res=cv::Mat_<float>::zeros(nrModes,m_chBasisVector.size()+1); // todo: change matrix shape internally, use STL vector?
	for (int mode=0; mode < nrModes; mode++) {
		nrXChannels = chCoeff.cols;
		xChannel = chCoeff.clone();
		for (int k=(int)m_chBasisVector.size()-1; k>=0; k--) {
			int nrYChannels = m_chBasisVector[k]->getNrChannels();
			nrXChannels /= nrYChannels;
			cv::Mat_<float> tmpCoeff(1, nrYChannels, 0.0f);
			// integrate all other dimensions
			for (int l=0; l<nrYChannels; l++)
				for (cv::MatIterator_<float> xCi=xChannel.begin()+l*nrXChannels; xCi!=xChannel.begin()+(l+1)*nrXChannels; ++xCi)
					tmpCoeff(l) += *xCi;
			cv::Mat_<float> pRes;
			if (k==(int)m_chBasisVector.size()-1) {
				(m_chBasisVector[k])->decode(pRes, tmpCoeff, mode+1);
				res(mode,k) = pRes(0,2*mode);
				(m_chBasisVector[k])->encode(tmpCoeff,pRes(cv::Range(0,1),cv::Range(2*mode,2*mode+1)));
				tmpCoeff *= pRes(0,2*mode+1);
			}
			else {
				m_chBasisVector[k]->decode(pRes, tmpCoeff, 1);
				res(mode,k) = pRes(0,0); // TODO: check whether second 0 must be 2*mode!!
				m_chBasisVector[k]->encode(tmpCoeff,pRes(cv::Range(0,1),cv::Range(0,1)));
				tmpCoeff *= pRes(0,1); // TODO: see two lines above
			}

			// backproject only mode "mode"
			int elem_counter = 0;
			for (cv::MatIterator_<float> tCi=tmpCoeff.begin(); tCi!=tmpCoeff.end(); ++tCi) {
				xChannel(cv::Range(0,1),cv::Range(nrXChannels*elem_counter,nrXChannels*(elem_counter+1))) *= *tCi;
				elem_counter++;
			}
			cv::sqrt(xChannel,xChannel);
			tmpCoeff = cv::Mat_<float>::zeros(1, nrXChannels);
			// integrate out latest dimension
			for (int r=0; r<nrYChannels; r++) {
				tmpCoeff += xChannel(cv::Range(0,1),cv::Range(r*nrXChannels, (r+1)*nrXChannels));
			}
			xChannel = tmpCoeff.mul(tmpCoeff,1.0/m_chBasisVector[0]->getNorm()); // xChannel should be the certainty in the end
		}
		res(mode,m_chBasisVector.size()) = xChannel(0)/m_chBasisVector[0]->getNorm();
	}
	res = res.reshape(0,1);
	return;
}


void cvl::CombinedChannelBasis::decode(cv::Mat_<float> &res, const cv::SparseMat_<float> &chCoeff, const int nrModes) const {
	int nrXChannels;
	cv::SparseMat_<float> xChannel;
	res=cv::Mat_<float>::zeros(nrModes,m_chBasisVector.size()+1); // todo: change matrix shape internally, use STL vector?
	for (int mode=0; mode < nrModes; mode++) {
		nrXChannels = chCoeff.size(0);
		xChannel = chCoeff.clone();
		for (int k=(int)m_chBasisVector.size()-1; k>=0; k--) {
			int nrYChannels = m_chBasisVector[k]->getNrChannels();
			nrXChannels /= nrYChannels;
			cv::Mat_<float> tmpCoeff(1, nrYChannels, 0.0f);
			// integrate all other dimensions
			for (cv::SparseMatIterator_<float> xCi=xChannel.begin(); xCi!=xChannel.end(); ++xCi)
					tmpCoeff((xCi.node()->idx[0])/nrXChannels) += *xCi;
			cv::Mat_<float> pRes;
			if (k==(int)m_chBasisVector.size()-1) {
				(m_chBasisVector[k])->decode(pRes, tmpCoeff, mode+1);
				res(mode,k) = pRes(0,2*mode);
				(m_chBasisVector[k])->encode(tmpCoeff,pRes(cv::Range(0,1),cv::Range(2*mode,2*mode+1)));
				tmpCoeff *= pRes(0,2*mode+1);
			}
			else {
				m_chBasisVector[k]->decode(pRes, tmpCoeff, 1);
				res(mode,k) = pRes(0,0); // TODO: check whether second 0 must be 2*mode!!
				m_chBasisVector[k]->encode(tmpCoeff,pRes(cv::Range(0,1),cv::Range(0,1)));
				tmpCoeff *= pRes(0,1); // TODO: see two lines above
			}

			// backproject only mode "mode"
			int size[] = {nrXChannels};
			cv::SparseMat_<float> tmpSCoeff(1,size);
			for (cv::SparseMatIterator_<float> xCi=xChannel.begin(); xCi!=xChannel.end(); ++xCi) {
				int YChIndex = (xCi.node()->idx[0])/nrXChannels;
				int XChIndex = (xCi.node()->idx[0])%nrXChannels;
				// integrate out latest dimension
				tmpSCoeff.ref(XChIndex) += sqrt(*xCi*tmpCoeff(YChIndex));
			}
			for (cv::SparseMatIterator_<float> tCi=tmpSCoeff.begin(); tCi!=tmpSCoeff.end(); ++tCi)
				*tCi *= *tCi/m_chBasisVector[0]->getNorm();  // xChannel should be the certainty in the end
			xChannel = tmpSCoeff;
		}
		res(mode,m_chBasisVector.size()) = xChannel(0)/m_chBasisVector[0]->getNorm();
	}
	res = res.reshape(0,1);
	return;
}

// ---------------   ChannelVector ---------------------

// Constructors

cvl::ChannelVector::ChannelVector(){m_chBasis=NULL;}

cvl::ChannelVector::ChannelVector(ChannelBasis *chBasis):cv::Mat_<float>(1,chBasis->getNrChannels(),0.0f){
	m_chBasis = chBasis;
	m_support.resize(2);
	m_support[0] = 1;
	m_support[1] = 1;
} 

cvl::ChannelVector::ChannelVector(ChannelBasis *chBasis, const cv::Mat_<float> coeffs):cv::Mat_<float>(coeffs){
	CV_Assert(chBasis->getNrChannels()==coeffs.cols);
		m_chBasis = chBasis;
		m_support.resize(2);
		m_support[0] = 1;
		m_support[1] = coeffs.rows;
}

void cvl::ChannelVector::setChannelBasis(ChannelBasis *chBasis){
	m_chBasis=chBasis;   
	create(1,m_chBasis->getNrChannels());
	setTo(0.0f);
	return;
}

void cvl::ChannelVector::addSample(const Mat &vals){
	// It is assumed that VALS is a k-channel 2D image (float), where k equals the number of channel dimensions
	// + [-2, -1, 0, 1]. The cases are as follows:
	// -2: CCFM is computed ... In the latter case, the last dimension indicates the certainty.
	// all elements must be float!
	cv::transpose(*this,*this); // todo: change matrix orientation internally!
	std::vector<int> sizes(1);
	m_chBasis->getNrChannelsVec(sizes);
	int nrK = vals.channels(), nrR = vals.rows, nrC = vals.cols, nrEls;
	int sizeDiff = nrK-sizes.size();
	bool createCCFM = false, dedicatedCertainty = false;
	switch (sizeDiff) {
	case -2:
		// two more channel dimensions means that we compute a CCFM without dedicated certainty
		createCCFM = true;
		dedicatedCertainty = false;
		break;
	case -1:
		// one more channel dimensions means that we compute a CCFM with dedicated certainty
		createCCFM = true;
		dedicatedCertainty = true;
		break;
	case 0:
		// same number of channels means that we compute pointwise channels without dedicated certainty
		createCCFM = false;
		dedicatedCertainty = false;
		break;
	case 1:
		// one less channels dimensions means that we compute pointwise channels with dedicated certainty
		createCCFM = false;
		dedicatedCertainty = true;
		break;
	default:
		// this must throw an error
		std::cerr << "Error: incompatible format of channel vector and data dimensionality! \n"; //FIXME: use CV_Error macros
		break;
	}
	//TODO: add check for continuous and loop more efficiently
	cv::Mat_<float> tmp_chCoeff;
	std::vector<float> cVals(sizes.size());
	if (createCCFM)
		nrEls = 1;
	else
		nrEls = nrR*nrC;
	if (rows != m_chBasis->getNrChannels() || cols != nrEls) {
		const int cSize[] = {m_chBasis->getNrChannels(), nrEls};
		create(2,cSize);
		setTo(0.0f);
		m_support.resize(2);
		if (createCCFM) {
			m_support[0] = 1;
			m_support[1] = 1;
		} else {
			m_support[0] = nrR;
			m_support[1] = nrC;
		}
	}
	for (int rdx = 0; rdx < nrR; rdx++) { //TODO: change to iterator
		const float* data = vals.ptr<float>(rdx);
		for (int cdx = 0; cdx < nrC; cdx++) {
			//TODO: this is slow, but I do not have another way to get the unknown number of elements into an STL-vector
			if (dedicatedCertainty)
				for (int kdx = 0; kdx < (nrK-1); kdx++)
					cVals[kdx] = data[cdx*nrK+kdx];
			else
				for (int kdx = 0; kdx < nrK; kdx++)
					cVals[kdx] = data[cdx*nrK+kdx];
			if (createCCFM) {
				cVals[sizes.size()-2] = rdx;
				cVals[sizes.size()-1] = cdx;
			}
			m_chBasis->encode(tmp_chCoeff,cVals);
			if (createCCFM){
				float corrW = pow(nrR*nrC,2.f/((float)sizes.size())); //TODO: check that normalization is correct
				if (dedicatedCertainty)
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						at<float>(kdx,0)+=tmp_chCoeff(kdx,0)*data[cdx*nrK+nrK-1]/corrW;
				else
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						at<float>(kdx,0)+=tmp_chCoeff(kdx,0)/corrW;
			}
			else
				if (dedicatedCertainty)
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						at<float>(kdx,rdx*nrC+cdx)+=tmp_chCoeff(kdx,0)*data[cdx*nrK+nrK-1];
				else
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						at<float>(kdx,rdx*nrC+cdx)+=tmp_chCoeff(kdx,0);
		}
	}
	cv::transpose(*this,*this);
	return;
}

void cvl::ChannelVector::normalize(){
	float nrm=norm(*this,cv::NORM_L1);
	if (nrm>0)
		(*this)*=(m_chBasis->getNorm()/nrm);
	else
		setTo(m_chBasis->getNorm()/cols);
	return;
}

void cvl::ChannelVector::channelImage(cv::Mat &res) const {
	res = reshape(cols,m_support[0]);
}

void cvl::ChannelVector::histogramMatrix(cv::Mat &res) const {
	std::vector<int> sizes(1);
	m_chBasis->getNrChannelsVec(sizes);
	sizes.push_back(rows);
	cv::Mat resTmp((int)sizes.size(),&sizes[0],CV_32F,data);
	res = resTmp;
}

void cvl::ChannelVector::decode(cv::Mat &res, const int nrModes) const {
	std::vector<int> sizes(1);
	m_chBasis->getNrChannelsVec(sizes);
	int nrChannels = nrModes*(sizes.size()+1);
	CV_Assert(nrChannels<=CV_CN_MAX);
	int nrEls = m_support[0]*m_support[1];
	cv::Mat_<float> localRes(nrEls,nrChannels);
	cv::Mat_<float> res_col(1,nrChannels);
	for (int cdx = 0; cdx<rows; cdx++) {
		m_chBasis->decode(res_col,row(cdx),nrModes);
		res_col.copyTo(localRes(cv::Range(cdx,cdx+1),cv::Range::all()));
	}
	res = localRes.reshape(nrChannels,m_support[0]);
}

std::ostream &operator << (std::ostream &os, const cvl::ChannelVector &v) {
	int nrR = v.rows;
	int nrC = v.cols;
	if (nrC>=0) {
		for (int j=0; j<nrR; j++) {
			os << std::endl << '[' << nrC << "](";
			if (nrC > 0)
				os << v.at<float>(j,0);
			for (int i = 1; i < nrC; ++ i)
				os << ',' << v.at<float>(j,i);
			os << ')';
		}
	}
	else
		os << "more than 2 dimensions";
	return os;
}

// ---------------   ChannelSVector ---------------------

// Constructors

cvl::ChannelSVector::ChannelSVector(){m_chBasis=NULL;}

cvl::ChannelSVector::ChannelSVector(ChannelBasis *chBasis):cv::SparseMat_<float>() {
	m_chBasis=chBasis;
	int sizes[] = {chBasis->getNrChannels()};
	create(1,sizes);
	m_support.resize(2);
	m_support[0] = 1;
	m_support[1] = 1;
} 

void cvl::ChannelSVector::setChannelBasis(ChannelBasis *chBasis){
	m_chBasis=chBasis;   
	int sizes[1] = {chBasis->getNrChannels()};
	create(1,sizes);
	return;
}

void cvl::ChannelSVector::addSample(const cv::Mat &vals){
	// It is assumed that VALS is a k-channel 2D image (float), where k equals the number of channel dimensions
	// + [-2, -1, 0, 1]. The cases are as follows:
	// -2: CCFM is computed ... In the latter case, the last dimension indicates the certainty.
	// all elements must be float!
	std::vector<int> sizes(1);
	m_chBasis->getNrChannelsVec(sizes);
	int nrK = vals.channels(), nrR = vals.rows, nrC = vals.cols, nrEls;
	int sizeDiff = nrK-sizes.size();
	bool createCCFM = false, dedicatedCertainty = false;
	switch (sizeDiff) {
	case -2:
		// two more channel dimensions means that we compute a CCFM without dedicated certainty
		createCCFM = true;
		dedicatedCertainty = false;
		break;
	case -1:
		// one more channel dimensions means that we compute a CCFM with dedicated certainty
		createCCFM = true;
		dedicatedCertainty = true;
		break;
	case 0:
		// same number of channels means that we compute pointwise channels without dedicated certainty
		createCCFM = false;
		dedicatedCertainty = false;
		break;
	case 1:
		// one less channels dimensions means that we compute pointwise channels with dedicated certainty
		createCCFM = false;
		dedicatedCertainty = true;
		break;
	default:
		// this must throw an error
		std::cerr << "Error: incompatible format of channel vector and data dimensionality! \n"; //FIXME: use CV_Error macros
		break;
	}
	//TODO: add check for continuous and loop more efficiently
	cv::SparseMat_<float> tmp_chCoeff;
	std::vector<float> cVals(sizes.size());
	if (createCCFM)
		nrEls = 1;
	else
		nrEls = nrR*nrC;
	if (dims() != 2 || size(0) != m_chBasis->getNrChannels() || size(1) != nrEls) {
		const int cSize[] = {m_chBasis->getNrChannels(), nrEls};
		create(2,cSize);
		clear();
		m_support.resize(2);
		if (createCCFM) {
			m_support[0] = 1;
			m_support[1] = 1;
		} else {
			m_support[0] = nrR;
			m_support[1] = nrC;
		}
	}
	for (int rdx = 0; rdx < nrR; rdx++) { //TODO: change to iterator
		const float* data = vals.ptr<float>(rdx);
		for (int cdx = 0; cdx < nrC; cdx++) {
			//TODO: this is slow, but I do not have another way to get the unknown number of elements into an STL-vector
			if (dedicatedCertainty)
				for (int kdx = 0; kdx < (nrK-1); kdx++)
					cVals[kdx] = data[cdx*nrK+kdx];
			else
				for (int kdx = 0; kdx < nrK; kdx++)
					cVals[kdx] = data[cdx*nrK+kdx];
			if (createCCFM) {
				cVals[sizes.size()-2] = rdx;
				cVals[sizes.size()-1] = cdx;
			}
			m_chBasis->encode(tmp_chCoeff,cVals);
			if (createCCFM){
				float corrW = pow(nrR*nrC,2.f/((float)sizes.size())); //TODO: check that normalization is correct
				if (dedicatedCertainty)
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						ref(kdx,0)+=tmp_chCoeff(kdx,0)*data[cdx*nrK+nrK-1]/corrW;
				else
					for (int kdx = 0; kdx < m_chBasis->getNrChannels(); kdx++)
						ref(kdx,0)+=tmp_chCoeff(kdx)/corrW;
			}
			else
				for (cv::SparseMatIterator_<float> tCi=tmp_chCoeff.begin(); tCi!=tmp_chCoeff.end(); ++tCi) {
					const Node* tCn = tCi.node();
					if (dedicatedCertainty)
						ref(tCn->idx[0],rdx*nrC+cdx) += *tCi*data[cdx*nrK+nrK-1];
					else
						ref(tCn->idx[0],rdx*nrC+cdx) += *tCi;
				}
		}
	}
	return;
}

void cvl::ChannelSVector::normalize(){
	float nrm=norm(*this,cv::NORM_L1);
	if (nrm>0) {
		for (cv::SparseMatIterator_<float> ti=begin(); ti!=end(); ++ti)
			*ti*=(m_chBasis->getNorm()/nrm);
	}
	else {
		int nchs = size(0);
		float val = m_chBasis->getNorm()/nchs;
		for (int n=0; n<nchs; ++n)
			ref(&n)=val;
	}
	return;
}

//Decoding
void cvl::ChannelSVector::decode(cv::Mat_<float> &res, const int nrModes) const {
	std::vector<int> sizes(1);
	m_chBasis->getNrChannelsVec(sizes);
	int nrChannels = nrModes*(sizes.size()+1);
	CV_Assert(nrChannels<=CV_CN_MAX);
	int nrEls = m_support[0]*m_support[1];
	res.create(nrEls,nrChannels);
	cv::Mat_<float> res_col(1,nrChannels);
//	cv::Mat tmpM;
//	convertTo(tmpM,CV_32F);
//	cv::Mat_<float> tmpM(size(1),size(0),0.0f); // todo: sparse vectors still have the wrong orientation
	std::vector<cv::SparseMat_<float> > tmpM(size(1));
	int sizeTmp[] = {size(0)};
	for (int cdx = 0; cdx<size(1); cdx++)
		tmpM[cdx].create(1,sizeTmp);
	for (cv::SparseMatConstIterator_<float> it = begin(); it != end(); ++it)
		(tmpM[it.node()->idx[1]]).ref(it.node()->idx[0])=*it;
	for (int cdx = 0; cdx<size(1); cdx++) {
		m_chBasis->decode(res_col,tmpM[cdx],nrModes);
		res_col.copyTo(res(cv::Range(cdx,cdx+1),cv::Range::all()));
	}
	res.reshape(nrChannels,m_support[0]);
}

std::ostream &operator << (std::ostream &os, const cvl::ChannelSVector &v) {
	int nrR = (int)v.size(0);
	os << '[' << nrR << "](";
	for (cv::SparseMatConstIterator_<float> i = v.begin(); i != v.end(); ++i) {
		const cv::SparseMat::Node* n = i.node();
		if (i == v.begin())
			os << n->idx[0] << ":" << *i;
		else
			os << ", " << n->idx[0] << ":" << *i;
	}
	os << ')';
	return os;
}
