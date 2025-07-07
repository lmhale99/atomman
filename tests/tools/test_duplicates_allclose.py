# coding: utf-8
from io import StringIO

# https://docs.pytest.org/en/latest/
import pytest

import atomman as am
import numpy as np
import pandas as pd

# csv of dataframe to use.
# first values are values to test against
# test_* values are the expected outputs for the tests
csv = """
str1,str2,int1,float1,float2,test_all_first,test_all_last,test_all_false,test_dcols_first,test_dcols_last,test_dcols_false,test_fcols_first,test_fcols_last,test_fcols_false
b,y,1,1.0448696892749498,10.043939753734117,False,False,False,False,True,True,False,False,False
a,y,1,0.9038123150145498,9.983496835984274,False,False,False,False,True,True,False,False,False
b,x,1,0.9803175864013501,10.014852172497893,False,False,False,False,True,True,False,False,False
b,x,0,0.9362356629957925,10.029225813896156,False,False,False,False,True,True,False,True,True
b,x,1,0.9520056087229036,10.049876177397062,False,False,False,True,True,True,False,False,False
b,x,0,1.0622617053862702,10.042542586852289,False,False,False,True,True,True,False,False,False
b,x,0,1.0869031817826926,10.038488740496314,False,False,False,True,True,True,False,False,False
a,x,0,1.01461464369603,9.989118972512804,False,False,False,False,True,True,False,False,False
b,y,0,1.0877286885396968,10.00178129242674,False,False,False,False,True,True,False,False,False
a,y,0,0.9066744165054729,9.951810227418571,False,False,False,False,True,True,True,False,True
b,y,0,0.9846141116833345,10.045992530965963,False,False,False,True,True,True,False,False,False
b,x,1,0.923652393892943,9.957550028291852,False,False,False,True,True,True,True,False,True
b,y,1,1.0196402820396742,10.029406071454774,False,True,True,True,True,True,False,False,False
a,y,0,1.004885491248553,10.003356909589773,False,False,False,True,True,True,False,False,False
a,y,0,0.9964793773923387,9.975273588507132,False,False,False,True,True,True,True,False,True
a,y,1,0.992184499467428,10.041090829702364,False,False,False,True,True,True,False,False,False
a,x,1,0.9215000472039374,9.952903941855064,False,False,False,False,True,True,False,True,True
a,x,0,0.9832887461209068,10.010110839683874,True,False,True,True,True,True,False,False,False
b,y,1,0.9406658316645825,10.040727850433978,False,False,False,True,True,True,True,False,True
b,y,1,0.9461438276327032,9.977564387629958,False,False,False,True,True,True,False,False,False
b,x,1,0.9387305946885973,10.039165302003397,False,False,False,True,True,True,True,True,True
b,y,1,1.0860973290856524,9.962448607332073,False,False,False,True,True,True,False,False,False
b,y,1,1.0073608049368545,10.044717350901736,False,False,False,True,True,True,False,False,False
b,y,1,0.9428629935232474,10.02417609441616,False,False,False,True,True,True,False,False,False
b,x,1,0.9814445645605034,9.9794411425119,False,False,False,True,True,True,False,False,False
b,y,1,1.029687084158268,10.033604760554907,True,False,True,True,True,True,False,True,True
a,x,0,1.0834798250665072,10.02556504335372,False,False,False,True,True,True,True,False,True
a,x,0,0.973661807593718,10.001197160305864,True,True,True,True,True,True,True,False,True
b,y,0,0.918981272623069,9.993458378057927,False,False,False,True,True,True,False,False,False
b,x,1,0.9891281633219332,9.950793952382199,False,False,False,True,False,True,False,False,False
a,x,1,0.9575685660206714,9.977286155494918,False,False,False,True,True,True,False,False,False
b,y,1,1.0541206289305616,9.986960917610544,False,False,False,True,True,True,False,False,False
b,y,1,1.0373252379911777,10.002411435039834,False,False,False,True,False,True,False,False,False
a,x,1,1.0032207995117206,9.98730533457245,False,False,False,True,True,True,False,False,False
a,x,0,0.9546262226748408,10.020232902905057,False,False,False,True,True,True,False,False,False
b,x,0,1.0267472790194225,9.970629801948604,False,False,False,True,True,True,False,False,False
a,y,1,1.0431515125318447,9.972156353460697,False,False,False,True,True,True,False,False,False
a,y,0,1.0725347870325252,10.026450441385988,False,False,False,True,True,True,False,True,True
b,x,0,0.9276425211456973,9.999500254100454,False,False,False,True,True,True,False,False,False
a,y,1,0.9044492486099063,9.951157892623481,False,False,False,True,True,True,False,True,True
b,y,0,0.9933505436366469,9.982308383731437,False,False,False,True,False,True,False,True,True
a,y,0,1.0331664038512356,10.03380086881164,False,False,False,True,True,True,True,False,True
a,x,1,0.9405401301453724,10.041333344863984,False,False,False,True,False,True,True,True,True
a,y,1,0.9814931412667157,9.955144807614696,False,False,False,True,False,True,False,False,False
b,x,0,0.903862052902848,10.042311950668246,False,False,False,True,True,True,False,False,False
a,x,0,0.9577648614340931,10.005398521486972,False,True,True,True,True,True,False,False,False
b,x,0,1.0622476123583193,9.969704710361052,False,False,False,True,False,True,False,False,False
a,y,0,0.9683425589137947,9.993520318883515,False,False,False,True,True,True,False,True,True
a,y,0,0.9822481430809309,10.036454435853742,False,False,False,True,False,True,False,False,False
a,x,0,1.058669936100845,10.025104728885823,False,False,False,True,False,True,False,False,False
"""

sliced_csv = """
,test_first,test_last,test_false
22,False,False,False
9,False,False,False
33,False,False,False
6,False,False,False
15,False,False,False
25,True,False,True
46,False,False,False
34,False,False,False
41,False,False,False
44,False,False,False
40,False,False,False
8,False,False,False
3,False,False,False
24,False,False,False
7,False,False,False
48,False,False,False
5,False,False,False
2,False,False,False
16,False,False,False
17,False,False,False
20,False,False,False
10,False,False,False
35,False,False,False
42,False,False,False
19,False,False,False
4,False,False,False
28,False,False,False
49,False,False,False
1,False,False,False
23,False,False,False
36,False,False,False
37,False,False,False
12,False,True,True
21,False,False,False
43,False,False,False
"""


class Test_duplicates_allclose():

    @property
    def dataframe(self):
        return pd.read_csv(StringIO(csv))
    
    @property
    def sliced_results(self):
        return pd.read_csv(StringIO(sliced_csv), index_col=0)

    def test_all(self):
        dcols = ['str1', 'str2', 'int1']
        fcols = {'float1':0.02, 'float2':0.01}
        
        df = self.dataframe

        test_all_first = am.tools.duplicated_allclose(df, dcols, fcols)
        assert np.all(test_all_first == df.test_all_first)

        test_all_first = am.tools.duplicated_allclose(df, dcols, fcols, keep='first')
        assert np.all(test_all_first == df.test_all_first)

        test_all_last = am.tools.duplicated_allclose(df, dcols, fcols, keep='last')
        assert np.all(test_all_last == df.test_all_last)

        test_all_false = am.tools.duplicated_allclose(df, dcols, fcols, keep=False)
        assert np.all(test_all_false == df.test_all_false)

        with pytest.raises(ValueError):
            am.tools.duplicated_allclose(df, dcols, fcols, keep='nope')

    def test_dcols(self):
        dcols = ['str1', 'str2', 'int1']
        fcols = {}
        
        df = self.dataframe

        test_dcols_first = am.tools.duplicated_allclose(df, dcols, fcols, keep='first')
        assert np.all(test_dcols_first == df.test_dcols_first)

        test_dcols_last = am.tools.duplicated_allclose(df, dcols, fcols, keep='last')
        assert np.all(test_dcols_last == df.test_dcols_last)

        test_dcols_false = am.tools.duplicated_allclose(df, dcols, fcols, keep=False)
        assert np.all(test_dcols_false == df.test_dcols_false)

    def test_fcols(self):
        dcols = []
        fcols = {'float1':0.02, 'float2':0.01}
        
        df = self.dataframe

        test_fcols_first = am.tools.duplicated_allclose(df, dcols, fcols, keep='first')
        assert np.all(test_fcols_first == df.test_fcols_first)

        test_fcols_last = am.tools.duplicated_allclose(df, dcols, fcols, keep='last')
        assert np.all(test_fcols_last == df.test_fcols_last)

        test_fcols_false = am.tools.duplicated_allclose(df, dcols, fcols, keep=False)
        assert np.all(test_fcols_false == df.test_fcols_false)

    def test_slice(self):
        dcols = ['str1', 'str2', 'int1']
        fcols = {'float1':0.02, 'float2':0.01}
        
        df = self.dataframe
        sliced_results = self.sliced_results

        # slice df according to indices of sliced results
        df2 = df.loc[sliced_results.index]

        test_first = am.tools.duplicated_allclose(df2, dcols, fcols, keep='first')
        assert np.all(test_first == sliced_results.test_first), test_first

        test_last = am.tools.duplicated_allclose(df2, dcols, fcols, keep='last')
        assert np.all(test_last == sliced_results.test_last)

        test_false = am.tools.duplicated_allclose(df2, dcols, fcols, keep=False)
        assert np.all(test_false == sliced_results.test_false)
