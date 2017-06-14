void SceECM::applyECMForce_M() {
        
        
	thrust::transform(nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount,
			  nodes->getInfoVecs().nodeLocY.begin()+ totalNodeCountForActiveCells,
			  nodes->getInfoVecs().nodeLocY.begin()) + allocPara_m.bdryNodeCount,
			AddECMForc(25.0);

}
