package filters;

import model.RoadModel;

import org.jblas.DoubleMatrix;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import doubleMatrix.Concat;
import doubleMatrix.InverseMatrix;
import section.*;


public abstract class Filter {

	public DoubleMatrix var;
	public DoubleMatrix mean;
	public DoubleMatrix csterm;
	public DoubleMatrix f_var;
	public DoubleMatrix f_mean;	
	public double gamma;
	public double gammaStar;
	public DoubleMatrix Kgain;	
	public DoubleMatrix measurements;
	public DoubleMatrix C1;
	public DoubleMatrix C2;
	public DoubleMatrix FG;
	public DoubleMatrix F;
	public DoubleMatrix G;	
	public DoubleMatrix Lambda;
	public DoubleMatrix Ltilde;
	public DoubleMatrix Lower;
	public DoubleMatrix I1;
	public DoubleMatrix I2;
	RoadModel roadModel;
	Estimation estimation;
	int size;
	int sizeMeasurements;
	int cellsec;	
	int stepMeasurements;
	DoubleMatrix modelVar;	
	DoubleMatrix measureVar;	
	DoubleMatrix priorVar;
	public DoubleMatrix scaling;
	double threshold;
	public DoubleMatrix measure; //public only for test	
	int numUp;	
	Section section;
	public double minLambda;
	public double maxLower;
	
	
	abstract public void nextStep();
	abstract public void nextStepNoData();
	
	protected void initial(Estimation _estimation) {
		
		estimation = _estimation;
		roadModel = estimation.roadModel;
		section =estimation.roadModel.section;
		getNewParametersFromModel();		
		mean = roadModel.initialMean;
		f_mean=mean.dup();
		var = roadModel.initialVar.mul(25);
		f_var=var.dup();
		priorVar=var.dup();		
		numUp = 0;
		Lambda=DoubleMatrix.zeros(size, size);
		csterm=DoubleMatrix.zeros(2, 1);
		scaling=DoubleMatrix.zeros(2, 1);
		scaling.put(0,0,1);
		scaling.put(1,0,1);
		threshold=0.01;
	}
	
	public static Filter createFilter(Estimation _estimation) {
		Filter f = null;
        f = new KCF(_estimation);
		return f;	 	
	}
	

	
	public void ComputeGammaStar(Filter _filter1, Filter _filter2){
		    double minbigLambda;
		    if (_filter1.minLambda<minLambda){
		    	minbigLambda=(_filter1.minLambda)/((double)3);
		    }
		    else{
		    	minbigLambda=(minLambda)/((double)3);
		    }
		    if (minbigLambda>(_filter2.minLambda)/((double)3)){
		    	minbigLambda=(_filter2.minLambda)/((double)3);
		    }
		    gammaStar=Math.sqrt(minbigLambda/maxLower);

	}
	public void ComputeGammaStar1(Filter _filter1){
	    double minbigLambda;
	    if ((_filter1.minLambda)<minLambda){
	    	minbigLambda=(_filter1.minLambda)/((double)3);
	    }
	    else{
	    	minbigLambda=(minLambda)/((double)3);
	    }
	    gammaStar=Math.sqrt(minbigLambda/maxLower);

}
	public void ComputeGammaStar2(Filter _filter2){
	    double minbigLambda;
	    if ((_filter2.minLambda)<minLambda){
	    	minbigLambda=(_filter2.minLambda)/((double)3);
	    }
	    else{
	    	minbigLambda=(minLambda)/((double)3);
	    }
	    gammaStar=Math.sqrt(minbigLambda/maxLower);

}
	
	public DoubleMatrix Consensus(Filter _filter1, Filter _filter2){
		int overlap=roadModel.trueSolution.overlapsec;
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) {
		    if(roadModel.section.getClass().getSimpleName().equals("FC")){
		    	return mean;
		    }
		    else{
		
				
				C1=new DoubleMatrix(size,overlap);
				C2=new DoubleMatrix(size,overlap);
				if (gammaStar>1){
					gammaStar=0.5;
				}
				
//				DoubleMatrix Xi1=DoubleMatrix.zeros(size,size);
//				DoubleMatrix Xi2=DoubleMatrix.zeros(size,size);
//				for (int i=0; i<overlap;i++){
//					Xi1.put(i,i,1);
//				}
//				for (int i=size-overlap;i<size;i++){
//					Xi2.put(i,i,1);
//				}			

				


//				C1=Xi1.mmul(var.mmul(I1).mul(gammaStar));
//				C2=Xi2.mmul(var.mmul(I2).mul(gammaStar));
				C1=var.mmul(I1).mul(gammaStar*0.99);
				C2=var.mmul(I2).mul(gammaStar*0.99);
				
				csterm.put(0,0,(C1.mmul(_filter1.mean.getRange(_filter1.size-overlap,_filter1.size,0,1).sub(mean.getRange(0,overlap,0,1)))).norm2());
				csterm.put(1,0,(C2.mmul(_filter2.mean.getRange(0,overlap,0,1).sub(mean.getRange(size-overlap,size,0,1)))).norm2());
				for (int i=0; i<2;i++){
					if (csterm.get(i)>threshold){
						scaling.put(i,0,((double)threshold)/((double)csterm.get(i)));
					}
					else{
						scaling.put(i,0,1);
					}
				}

				
				DoubleMatrix mean1=new DoubleMatrix(size,1);
				mean1=mean.add((C1.mmul(_filter1.mean.getRange(_filter1.size-overlap,_filter1.size,0,1).sub(mean.getRange(0,overlap,0,1)))).mmul(scaling.get(0))).add((C2.mmul(_filter2.mean.getRange(0,overlap,0,1).sub(mean.getRange(size-overlap,size,0,1)))).mmul(scaling.get(1)));

				
				return mean1;
				
		
		    }

			}
			else{return mean;}
	}
	
    public DoubleMatrix Consensus1(Filter _filter1){
    	int overlap=roadModel.trueSolution.overlapsec;
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) {
		    if(roadModel.section.getClass().getSimpleName().equals("FC")){
		    	return mean;
		    }
		    else{
		    	
				C1=new DoubleMatrix(size,overlap);
				if (gammaStar>1){
					gammaStar=0.5;
				}
				
//				DoubleMatrix Xi1=DoubleMatrix.zeros(size,size);
//				DoubleMatrix Xi2=DoubleMatrix.zeros(size,size);
//				for (int i=0; i<overlap;i++){
//					Xi1.put(i,i,1);
//				}
//				for (int i=size-overlap;i<size;i++){
//					Xi2.put(i,i,1);
//				}			



//				C1=Xi1.mmul(var.mmul(I1).mul(gammaStar));
				C1=var.mmul(I1).mul(gammaStar*0.99);
				
				csterm.put(0,0,(C1.mmul(_filter1.mean.getRange(_filter1.size-overlap,_filter1.size,0,1).sub(mean.getRange(0,overlap,0,1)))).norm2());
				
				for (int i=0; i<1;i++){
					if (csterm.get(i)>threshold){
						scaling.put(i,0,((double)threshold)/((double)csterm.get(i)));
					}
					else{
						scaling.put(i,0,1);
					}
				}
			
				DoubleMatrix mean1=new DoubleMatrix(size,1);
				mean1=mean.add((C1.mmul(_filter1.mean.getRange(_filter1.size-overlap,_filter1.size,0,1).sub(mean.getRange(0,overlap,0,1)))).mmul(csterm.get(0)));

				return mean1;
		    }

		}
		else {return mean;}
	}
	
    public DoubleMatrix Consensus2(Filter _filter2){
    	int overlap=roadModel.trueSolution.overlapsec;
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) {
		    if(roadModel.section.getClass().getSimpleName().equals("FC")){
		    	return mean;
		    }
		    else{
	
				C2=new DoubleMatrix(size,overlap);
				
				if (gammaStar>1){
					gammaStar=0.5;
				}
//				DoubleMatrix Xi1=DoubleMatrix.zeros(size,size);
//				DoubleMatrix Xi2=DoubleMatrix.zeros(size,size);
//				for (int i=0; i<overlap;i++){
//					Xi1.put(i,i,1);
//				}
//				for (int i=size-overlap;i<size;i++){
//					Xi2.put(i,i,1);
//				}			

				
//				C2=Xi2.mmul(var.mmul(I2).mul(gammaStar));
				C2=var.mmul(I2).mul(gammaStar*0.99);
				
				csterm.put(1,0,(C2.mmul(_filter2.mean.getRange(_filter2.size-overlap,_filter2.size,0,1).sub(mean.getRange(0,overlap,0,1)))).norm2());			
				for (int i=1; i<2;i++){
					if (csterm.get(i)>threshold){
						scaling.put(i,0,((double)threshold)/((double)csterm.get(i)));
					}
					else{
						scaling.put(i,0,1);
					}
				}
				
				DoubleMatrix mean1=new DoubleMatrix(size,1);
				mean1=mean.add((C2.mmul(_filter2.mean.getRange(_filter2.size-overlap,_filter2.size,0,1).sub(mean.getRange(0,overlap,0,1)))).mmul(csterm.get(1)));

				return mean1;
		    }

		}
		else {return mean;}
	}
    
	protected DoubleMatrix propagate(DoubleMatrix _density) {
		return estimation.propagate(_density);
	}	
	protected DoubleMatrix getMeasurements() {
		return roadModel.getMeasureVector();
	}	
	protected DoubleMatrix computeVar(DoubleMatrix _Var){
		DoubleMatrix _fVar=DoubleMatrix.zeros(_Var.getRows(), _Var.getColumns());
		_fVar=roadModel.section.ModelA.mmul(_Var).mmul(roadModel.section.ModelA.transpose()).add(modelVar);
		return _fVar;
	}	
	protected void getNewParametersFromModel() {
		size = roadModel.size;
		modelVar = roadModel.modelVar;
		sizeMeasurements = roadModel.sizeMeasurements;
		stepMeasurements = roadModel.stepMeasurements;
		measureVar = roadModel.measureVar;
		measure = roadModel.measure;
	}
	public void ComputeEig(){
		
				DoubleMatrix I=DoubleMatrix.eye(size);
				DoubleMatrix S=measure.transpose().mmul(InverseMatrix.invPoSym(measureVar)).mmul(measure);
				F=I.sub(var.mmul(S));
//				F=I;
				double mu=1;

				Lambda=InverseMatrix.invPoSym(F.mmul(roadModel.section.ModelA).mmul(priorVar.mul(mu*mu)).mmul(roadModel.section.ModelA.transpose()).mmul(F.transpose())).sub(InverseMatrix.invPoSym((F.mmul(roadModel.section.ModelA).mmul(priorVar.mul(mu*mu)).mmul(roadModel.section.ModelA.transpose()).mmul(F.transpose())).add(F.mmul((f_var.mmul(S.mul(mu*mu)).mmul(f_var)).add(roadModel.modelVar.mul(mu*mu))).mmul(F.transpose()))));
//				Lambda=InverseMatrix.invPoSym((roadModel.section.ModelA).mmul(priorVar).mmul(roadModel.section.ModelA.transpose())).sub(InverseMatrix.invPoSym(((roadModel.section.ModelA).mmul(priorVar).mmul(roadModel.section.ModelA.transpose())).add((f_var.mmul(S).mmul(f_var)))));	

				
				
				if (roadModel.trueSolution.index==0){
					
					Lower=Ltilde.transpose().mmul(I2.transpose()).mmul(var.mmul(I2)).mmul(Ltilde);
				}
				else if (roadModel.trueSolution.index==roadModel.trueSolution.numSections-1){
					Lower=Ltilde.transpose().mmul(I1.transpose()).mmul(var.mmul(I1)).mmul(Ltilde);
				}
				else{
					Lower=Ltilde.transpose().mmul((DoubleMatrix.concatHorizontally(I1, I2)).transpose()).mmul(var.mmul((DoubleMatrix.concatHorizontally(I1, I2)))).mmul(Ltilde);
				}
				double[][] _Lambda=new double[Lambda.rows][Lambda.columns];
				double[][] _Lower=new double[Lower.rows][Lower.columns];
				for (int i=0;i<Lambda.rows;i++){
					for (int j=0;j<Lambda.columns;j++){
						_Lambda[i][j]=Lambda.get(i, j);
					}
				}
				for (int i=0;i<Lower.rows;i++){
					for (int j=0;j<Lower.columns;j++){
						_Lower[i][j]=Lower.get(i, j);
					}
				}
				Matrix MLambda=new Matrix(_Lambda);
				Matrix MLower=new Matrix(_Lower);
				EigenvalueDecomposition MLambdaEig=new EigenvalueDecomposition(MLambda);
				EigenvalueDecomposition MLowerEig=new EigenvalueDecomposition(MLower);
				Matrix _MLambdaEig=MLambdaEig.getD();
				Matrix _MLowerEig=MLowerEig.getD();					
				minLambda=_MLambdaEig.get(0, 0);
				for(int i=1;i<_MLambdaEig.getColumnDimension();i++){
					if(_MLambdaEig.get(i, i)<minLambda){
						minLambda=_MLambdaEig.get(i, i);
					}						
				}					
				maxLower=_MLowerEig.get(0, 0);
				for(int i=1;i<_MLowerEig.getColumnDimension();i++){
					if(_MLowerEig.get(i, i)>maxLower){
						maxLower=_MLowerEig.get(i, i);
					}							
				}
		 }
					

	
}
