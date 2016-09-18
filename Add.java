package FPU;

/**
 * Created by mohammad on 6/28/2016.
 */
public class Add {

    String returnVal = "1";
    int bits = 23, bias = 127, Total = 32;
    Structure Alpha = new Structure();
    Structure Beta  = new Structure();
    Structure Gamma = new Structure();
    Structure expectedStruct = new Structure();
    public Structure IEEERepresentation(Structure structure, float TheNum, int bits, int bias){
        //	if (TheNum==0) {
        structure.IEEE=0;
        structure.error=0;
        structure.maxerror=0;
        structure.Sign= false;
        structure.realSign=0;
        structure.f=0;
        structure.e=0;
        structure.ZeroF=1;
        structure.InfF=0;
        structure.DenF=0;
        structure.NanF= 0;
        structure.OveF= false;
        structure.UndF=false;
        structure.G=0;
        structure.R=0;
        structure.S=0;


        if (TheNum!=0) {
            structure.Sign = (TheNum < 0);
            structure.realSign = (TheNum < 0) ? -1 : 1;
            TheNum=Math.abs(TheNum);

            float Y= (float) (Math.log(TheNum)/Math.log(2));
            float N= (float) Math.floor(Y);
            float Z=Y-N;

            float ffull= (float) (TheNum*Math.pow(2,-N));

            if (ffull >= 2){
                ffull = ffull/2;
                N = N + 1;
            }
            ffull = ffull - 1;
            float fdecimal= (float) Math.floor(ffull*Math.pow(2,bits));

            structure.f= (float) (fdecimal*Math.pow(2,-bits));
            structure.e= (int) (bias+N);
            structure.IEEE= (float) (structure.realSign * (1+structure.f)*Math.pow(2,structure.e-bias));
            structure.UndF = false;
            structure.DenF = 0;


            if (TheNum < Math.pow(2,-bias+1)){
                N = -bias+1;
                ffull= (float) (TheNum/Math.pow(2,N));
                fdecimal= (float) Math.floor(ffull*Math.pow(2,bits));
                structure.f= (float) (fdecimal/Math.pow(2,bits));
                structure.e= 1;
                structure.DenF = 1;
                structure.IEEE= (float) (structure.realSign*structure.f*Math.pow(2,N));
            }


            if (TheNum < Math.pow(2,-bias+1-bits)){
                structure.DenF = 0;
                structure.UndF = false;
            }

            structure.error=TheNum-Math.abs(structure.IEEE);
            structure.maxerror= (float) Math.pow(2,structure.e-bias-bits);


            structure.ZeroF=0;
            structure.InfF = 0;
            structure.NanF = 0;
            if (structure.e==2*bias+1){
                if (structure.f==0){
                    structure.InfF = 1;
                } else {
                    structure.NanF = 1;
                }
            }
            structure.OveF = (structure.e > 2*bias+1);
        }
//        structure.e_str = unsignedToBinStr(Structure.e, Total-bits-1 );
//        structure.f_str = unsignedToBinStr(Structure.bits);
        return structure;
    }
//
//    public String unsignedToBinStr(float TheNum, int bits) {
//        // initialize the string
//        String S = "";
//        float CurPower;
//        for (int CurPlace = bits-1; CurPlace >= 0; CurPlace--) {
//            CurPower = (float) Math.pow(2, CurPlace);
//            if (TheNum >= CurPower) {
//                S += "1";
//                TheNum = TheNum - CurPower;
//            } else  {
//                S += "0";
//            }
//        }
//        return S;
//    }
//    public void  signedToBinStr(float TheNum, int bits) {
//        int mark = 0;
//        if (TheNum < 0) {
//            mark = 1;
//            TheNum = -1 * TheNum;
//        }
//        if (mark == 1) {
//            // convert to 2's complements
//            for (int i = bits-1; i>=0; i--) {
//                if (state == 0) { // before first 1
//                    Res = S.charAt(i) + Res;
//                    if (S.charAt(i) == '1')
//                        state = 1;
//                }
//                else {
//                    Res = (S.charAt(i)=='1'?'0':'1') + Res;
//                }
//            }
//        }
//    }


    private float Add(Structure Struct1, Structure Struct2, Structure Struct3, Structure Struct4){
        // Is A = 0 ?
        float f1;
        if (Struct1.ZeroF != 0 || Struct1.DenF != 0){
            f1 = Struct1.f;
        }else{
            f1 = Struct1.f + 1;
        }
        // Is B = 0 ?
        float f2;
        if (Struct2.ZeroF != 0 || Struct2.DenF != 0){
            f2 = Struct2.f;
        }else{
            f2 = Struct2.f + 1;
        }

        int e1 = Struct1.e;
        int e2 = Struct2.e;

        /***************alignment step***********/
        int bitsDiff = 0;
        if (e1 > e2) {
            bitsDiff = (e1 - e2);
            f2 = (float) (f2 / Math.pow(2,bitsDiff));
            e2 = e1;
        } else {
            if (e1 < e2) {
                bitsDiff = (e2 - e1);
                f1 = (float) (f1 / Math.pow(2,bitsDiff));
                e1 = e2;
            } else {
                bitsDiff=0;
            }
        }

        float Fullf1 = f1;
        float Fullf2 = f2;

        f1 = TruncateToBits(f1,bits);
        f2 = TruncateToBits(f2,bits);
        float Extra1 = (float) ((Fullf1 - f1)*Math.pow(2,bits));
        float Extra2 = (float) ((Fullf2 - f2)*Math.pow(2,bits));

        Struct1.G = 0;
        if (Extra1 >=0.5 ){
            Struct1.G = 1;
            Extra1 = (float) (Extra1 - 0.5);
        }

        Struct2.G = 0;
        if (Extra2 >=0.5 ){
            Struct2.G = 1;
            Extra2 = (float) (Extra2 - 0.5);
        }

        Struct1.R = 0;
        if (Extra1 >= 0.25 ){
            Struct1.R = 1;
            Extra1 = (float) (Extra1 - 0.25);
        }

        Struct2.R = 0;
        if (Extra2 >= 0.25){
            Struct2.R = 1;
            Extra2 = (float) (Extra2 - 0.25);
        }

        Struct1.S = 0;
        if (Extra1 > 0 ){
            Struct1.S = 1;
        }

        Struct2.S = 0;
        if (Extra2 > 0 ){
            Struct2.S = 1;
        }

        int carryin = 0;

        int f1sign = Struct1.realSign;
        int f2sign = Struct2.realSign;

        float mna1= f1sign*f1;

        float mna2= f2sign*f2;

        /*TODO: check the != 0 condition*/
        if(!(Struct1.G != 0|| Struct1.R != 0 || Struct1.S != 0) && !(Struct2.G != 0|| Struct2.R !=0|| Struct2.S != 0))
        {
        }
        else
        {
            if(((f1sign != f2sign) && (f1 > f2) && !(Struct1.G != 0|| Struct1.R != 0|| Struct1.S != 0)))
            {
                carryin = (int) Math.pow(2,-23);

                if(f1 < f2)
                {
                    if(f2sign == -1)
                        mna2 = f2sign*f2 + carryin;
                    else
                        mna2 = f2sign*f2 - carryin;

                }
                else
                {
                    if(f1sign == -1)
                        mna1 = f1sign*f1 + carryin;
                    else
                        mna1 = f1sign*f1 - carryin;
                }
            }
            else if(((f2sign != f1sign)&& (f2 > f1) && !(Struct2.G != 0|| Struct2.R != 0|| Struct2.S != 0)))
            {
                carryin = (int) Math.pow(2,-23);
                if(f1 < f2)
                {
                    if(f2sign == -1)
                        mna2 = f2sign*f2 + carryin;
                    else
                        mna2 = f2sign*f2 - carryin;
                }
                else
                {
                    if(f1sign == -1)
                        mna1 = f1sign*f1 + carryin;
                    else
                        mna1 = f1sign*f1 - carryin;
                }

            }
        }

        float f3 = mna1 + mna2;
        float TheNum= (float) (f3*Math.pow(2,e1-bias));
        Struct3 = IEEERepresentation(Struct3, TheNum, bits, bias);

        boolean f3sign = (f3 < 0);
        int f3realSign = (f3 < 0) ? -1 : 1;
        f3 = Math.abs(f3);

        // Find G, R ans S for f3
        float Fullf3 = f3;
        f3 = TruncateToBits(f3,bits);
        float Extra3 = (float) ((Fullf3 - f3)*Math.pow(2,bits));

        //f3orig = Struct1.realSign*f1 + Math.pow(-1,Struct2.Sign^Subtract)*f2;
        int asign = 0;
        int bsign = 0;
        // post calculation modification
//	    if(Subtract)
        //          GRSf3=Struct1.realSign*(Struct1.G*4 + Struct1.R*2 + Struct1.S)- Struct2.realSign*(Struct2.G*4+Struct2.R*2+Struct2.S);
        //    else

        if(!(Struct1.G != 0|| Struct1.R != 0|| Struct1.S != 0) && !(Struct2.G != 0|| Struct2.R != 0|| Struct2.S !=0))
        {
        }
        else
        {

            if(((f1sign != f2sign ) && (f1 > f2) && !(Struct1.G != 0|| Struct1.R != 0|| Struct1.S !=0)))
            {
                asign = 1;
            }
            else if(((f2sign != f1sign)&& (f2 > f1) && !(Struct2.G != 0|| Struct2.R != 0|| Struct2.S != 0)))
            {
                bsign = 1;
            }

        }

        int GRSf3= (int) (f1sign*(asign*8 + Struct1.G*4 + Struct1.R*2 + Struct1.S)+ f2sign*(bsign*8 + Struct2.G*4+Struct2.R*2+Struct2.S));

        if (GRSf3<0)
        {
            GRSf3=Math.abs(GRSf3) & 7;
        }
        int Gf3=0;
        int Rf3=0;
        int Sf3=0;
        if (GRSf3 >= 4)
        {
            Gf3=1;
            GRSf3-=4;
        }
        if (GRSf3 >=2)
        {
            Rf3=1;
            GRSf3-=2;
        }
        if (GRSf3 == 1)
        {
            Sf3=1;
        }
        f3 = (float) (f3 + Gf3*Math.pow(2,-bits-1) +  Rf3*Math.pow(2,-bits-2)
                        + Sf3*Math.pow(2,-bits-3));
        float result = (float) (f3realSign*f3*Math.pow(2,(e2 - bias)));
        Struct3 = IEEERepresentation(Struct3, result, bits+3, bias);

        Fullf3 = Struct3.f;
        float e3 = Struct3.e;
        float R1, R2;
        if(Math.abs(f3)<1 && e1!=e2)
        {
            f3+=Gf3*Math.pow(2,-bias);
            R1=Rf3;
            R2=Sf3;
        }else{
            R1=Gf3;
            R2=(Sf3 != 0 || Rf3 != 0)? 1:0;
        }
        Struct3 = IEEERepresentation(Struct3, result, bits, bias);
        f3 = Struct3.f;

        if (Struct3.DenF != 0.0f || Struct3.ZeroF != 0){
            f3 = f3;
        } else {
            f3 = f3 + 1;
        }

        // ************************Find Final result using rounding techniques and f3, R1, R2
        // Truncation
        float f = f3;
        float e = e3;

        result = (float) (f3realSign*f*Math.pow(2,e-bias));
        Struct3 = IEEERepresentation(Struct3, result, bits, bias);


        // Round to nearest Even
        float tmp = 0;
        if ( R1 != 0 && R2 != 0) {
            tmp = R1;
        }
        else if (R1 == 0 || R2 == 0) {
            tmp = 0;
        }
        /*TODO: replaced with another calue check it out*/
        f = (float) (f3 + (tmp)*Math.pow(2,-bits));
        float temp = (float) ((f3*Math.pow(2,bits-1)) - Math.floor(f3*Math.pow(2,bits-1)));
        int lastBit = 0;
        if (temp >= 0.5  ) { lastBit = 1;}

        float feven = f;
//	fodd = f;
        if ((R1 == 1) && (R2 == 0)) { feven = (float) (feven + lastBit*Math.pow(2,-bits));
        }
//	if ((R1 == 1) && (R2 == 0)) { fodd = fodd + (lastBit^1)*pow(2,-bits)}
        e = e3;
        result = (float) (f3realSign*feven*Math.pow(2,e-bias));
        Struct3 = IEEERepresentation(Struct3, result, bits, bias);

        Struct4.f = Struct3.f;
        Struct4.e = Struct3.e;
        Struct4.IEEE = Struct3.IEEE;

        float MyResult = Struct3.IEEE;

        return MyResult;

    }
    public float  TruncateToBits(float Number, int bits){
        if (Number >= 0) {
            return (float) (Math.floor(Number * Math.pow(2,bits)) / Math.pow(2,bits));
        } else {
            return (float) (-Math.floor(-Number * Math.pow(2,bits)) / Math.pow(2,bits));
        }
    }
    public void FPAdder(float AVal, float BVal){
        if (Float.isNaN(AVal)){
            //return error massage
            System.out.println("A is an invalid number");
            return;
        }
        if (Float.isNaN(BVal)) {
            //return error massage
            System.out.println("A is an invalid number");
            return;
        }

        Alpha = IEEERepresentation(Alpha, AVal, bits, bias);
        Beta  = IEEERepresentation(Beta, BVal, bits, bias);

        // Before Add/Sub Check for overflow, underflow, zero, Inf, NaN
        if (Alpha.OveF){
            System.out.println(" An Overflow Occured in A");
            return;
        }
        if (Alpha.UndF) {
            System.out.println(" An underflow Occured in A");
            return;
        }
        if (Beta.OveF){
            System.out.println(" An Overflow Occured in B");
            return;
        }
        if (Beta.UndF) {
            System.out.println(" An UnderFlow Occured in B");
            return;
        }


        Structure returnStruct = new Structure();
        float res = Add(Alpha, Beta, Gamma,returnStruct);

        System.out.println(res);
    }

}


class Structure{

    public float IEEE;
    public float error;
    public float maxerror;
    public boolean Sign;
    public float f;
    public String f_str;
    public int e;
    public String e_str;
    public int ZeroF;
    public float InfF;
    public float DenF;
    public float NanF;
    public boolean OveF;
    public boolean UndF;
    public float G;
    public float R;
    public float S;
    public int realSign;

}