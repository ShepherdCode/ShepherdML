Уа+
бЃ
8
Const
output"dtype"
valuetensor"
dtypetype

NoOp
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype
О
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring 

VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 "serve*2.3.02v2.3.0-rc2-23-gb36436b0878еЧ)

embedding_2/embeddingsVarHandleOp*
_output_shapes
: *
dtype0*
shape:	*'
shared_nameembedding_2/embeddings

*embedding_2/embeddings/Read/ReadVariableOpReadVariableOpembedding_2/embeddings*
_output_shapes
:	*
dtype0
x
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*
shared_namedense_2/kernel
q
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes

:*
dtype0
p
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_2/bias
i
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes
:*
dtype0
f
	Adam/iterVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	Adam/iter
_
Adam/iter/Read/ReadVariableOpReadVariableOp	Adam/iter*
_output_shapes
: *
dtype0	
j
Adam/beta_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_1
c
Adam/beta_1/Read/ReadVariableOpReadVariableOpAdam/beta_1*
_output_shapes
: *
dtype0
j
Adam/beta_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/beta_2
c
Adam/beta_2/Read/ReadVariableOpReadVariableOpAdam/beta_2*
_output_shapes
: *
dtype0
h

Adam/decayVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name
Adam/decay
a
Adam/decay/Read/ReadVariableOpReadVariableOp
Adam/decay*
_output_shapes
: *
dtype0
x
Adam/learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/learning_rate
q
&Adam/learning_rate/Read/ReadVariableOpReadVariableOpAdam/learning_rate*
_output_shapes
: *
dtype0

gru_4/gru_cell_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*(
shared_namegru_4/gru_cell_4/kernel

+gru_4/gru_cell_4/kernel/Read/ReadVariableOpReadVariableOpgru_4/gru_cell_4/kernel*
_output_shapes

:0*
dtype0

!gru_4/gru_cell_4/recurrent_kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*2
shared_name#!gru_4/gru_cell_4/recurrent_kernel

5gru_4/gru_cell_4/recurrent_kernel/Read/ReadVariableOpReadVariableOp!gru_4/gru_cell_4/recurrent_kernel*
_output_shapes

:0*
dtype0

gru_4/gru_cell_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*&
shared_namegru_4/gru_cell_4/bias

)gru_4/gru_cell_4/bias/Read/ReadVariableOpReadVariableOpgru_4/gru_cell_4/bias*
_output_shapes

:0*
dtype0

gru_5/gru_cell_5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*(
shared_namegru_5/gru_cell_5/kernel

+gru_5/gru_cell_5/kernel/Read/ReadVariableOpReadVariableOpgru_5/gru_cell_5/kernel*
_output_shapes

:0*
dtype0

!gru_5/gru_cell_5/recurrent_kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*2
shared_name#!gru_5/gru_cell_5/recurrent_kernel

5gru_5/gru_cell_5/recurrent_kernel/Read/ReadVariableOpReadVariableOp!gru_5/gru_cell_5/recurrent_kernel*
_output_shapes

:0*
dtype0

gru_5/gru_cell_5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*&
shared_namegru_5/gru_cell_5/bias

)gru_5/gru_cell_5/bias/Read/ReadVariableOpReadVariableOpgru_5/gru_cell_5/bias*
_output_shapes

:0*
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0

Adam/embedding_2/embeddings/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:	*.
shared_nameAdam/embedding_2/embeddings/m

1Adam/embedding_2/embeddings/m/Read/ReadVariableOpReadVariableOpAdam/embedding_2/embeddings/m*
_output_shapes
:	*
dtype0

Adam/dense_2/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_2/kernel/m

)Adam/dense_2/kernel/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/m*
_output_shapes

:*
dtype0
~
Adam/dense_2/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_2/bias/m
w
'Adam/dense_2/bias/m/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/m*
_output_shapes
:*
dtype0

Adam/gru_4/gru_cell_4/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*/
shared_name Adam/gru_4/gru_cell_4/kernel/m

2Adam/gru_4/gru_cell_4/kernel/m/Read/ReadVariableOpReadVariableOpAdam/gru_4/gru_cell_4/kernel/m*
_output_shapes

:0*
dtype0
Ќ
(Adam/gru_4/gru_cell_4/recurrent_kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*9
shared_name*(Adam/gru_4/gru_cell_4/recurrent_kernel/m
Ѕ
<Adam/gru_4/gru_cell_4/recurrent_kernel/m/Read/ReadVariableOpReadVariableOp(Adam/gru_4/gru_cell_4/recurrent_kernel/m*
_output_shapes

:0*
dtype0

Adam/gru_4/gru_cell_4/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*-
shared_nameAdam/gru_4/gru_cell_4/bias/m

0Adam/gru_4/gru_cell_4/bias/m/Read/ReadVariableOpReadVariableOpAdam/gru_4/gru_cell_4/bias/m*
_output_shapes

:0*
dtype0

Adam/gru_5/gru_cell_5/kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*/
shared_name Adam/gru_5/gru_cell_5/kernel/m

2Adam/gru_5/gru_cell_5/kernel/m/Read/ReadVariableOpReadVariableOpAdam/gru_5/gru_cell_5/kernel/m*
_output_shapes

:0*
dtype0
Ќ
(Adam/gru_5/gru_cell_5/recurrent_kernel/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*9
shared_name*(Adam/gru_5/gru_cell_5/recurrent_kernel/m
Ѕ
<Adam/gru_5/gru_cell_5/recurrent_kernel/m/Read/ReadVariableOpReadVariableOp(Adam/gru_5/gru_cell_5/recurrent_kernel/m*
_output_shapes

:0*
dtype0

Adam/gru_5/gru_cell_5/bias/mVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*-
shared_nameAdam/gru_5/gru_cell_5/bias/m

0Adam/gru_5/gru_cell_5/bias/m/Read/ReadVariableOpReadVariableOpAdam/gru_5/gru_cell_5/bias/m*
_output_shapes

:0*
dtype0

Adam/embedding_2/embeddings/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:	*.
shared_nameAdam/embedding_2/embeddings/v

1Adam/embedding_2/embeddings/v/Read/ReadVariableOpReadVariableOpAdam/embedding_2/embeddings/v*
_output_shapes
:	*
dtype0

Adam/dense_2/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:*&
shared_nameAdam/dense_2/kernel/v

)Adam/dense_2/kernel/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/kernel/v*
_output_shapes

:*
dtype0
~
Adam/dense_2/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape:*$
shared_nameAdam/dense_2/bias/v
w
'Adam/dense_2/bias/v/Read/ReadVariableOpReadVariableOpAdam/dense_2/bias/v*
_output_shapes
:*
dtype0

Adam/gru_4/gru_cell_4/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*/
shared_name Adam/gru_4/gru_cell_4/kernel/v

2Adam/gru_4/gru_cell_4/kernel/v/Read/ReadVariableOpReadVariableOpAdam/gru_4/gru_cell_4/kernel/v*
_output_shapes

:0*
dtype0
Ќ
(Adam/gru_4/gru_cell_4/recurrent_kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*9
shared_name*(Adam/gru_4/gru_cell_4/recurrent_kernel/v
Ѕ
<Adam/gru_4/gru_cell_4/recurrent_kernel/v/Read/ReadVariableOpReadVariableOp(Adam/gru_4/gru_cell_4/recurrent_kernel/v*
_output_shapes

:0*
dtype0

Adam/gru_4/gru_cell_4/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*-
shared_nameAdam/gru_4/gru_cell_4/bias/v

0Adam/gru_4/gru_cell_4/bias/v/Read/ReadVariableOpReadVariableOpAdam/gru_4/gru_cell_4/bias/v*
_output_shapes

:0*
dtype0

Adam/gru_5/gru_cell_5/kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*/
shared_name Adam/gru_5/gru_cell_5/kernel/v

2Adam/gru_5/gru_cell_5/kernel/v/Read/ReadVariableOpReadVariableOpAdam/gru_5/gru_cell_5/kernel/v*
_output_shapes

:0*
dtype0
Ќ
(Adam/gru_5/gru_cell_5/recurrent_kernel/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*9
shared_name*(Adam/gru_5/gru_cell_5/recurrent_kernel/v
Ѕ
<Adam/gru_5/gru_cell_5/recurrent_kernel/v/Read/ReadVariableOpReadVariableOp(Adam/gru_5/gru_cell_5/recurrent_kernel/v*
_output_shapes

:0*
dtype0

Adam/gru_5/gru_cell_5/bias/vVarHandleOp*
_output_shapes
: *
dtype0*
shape
:0*-
shared_nameAdam/gru_5/gru_cell_5/bias/v

0Adam/gru_5/gru_cell_5/bias/v/Read/ReadVariableOpReadVariableOpAdam/gru_5/gru_cell_5/bias/v*
_output_shapes

:0*
dtype0

NoOpNoOp
А8
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*ы7
valueс7Bо7 Bз7

layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
	optimizer
trainable_variables
regularization_losses
	variables
		keras_api


signatures
b

embeddings
trainable_variables
regularization_losses
	variables
	keras_api
l
cell

state_spec
trainable_variables
regularization_losses
	variables
	keras_api
l
cell

state_spec
trainable_variables
regularization_losses
	variables
	keras_api
h

kernel
bias
trainable_variables
regularization_losses
 	variables
!	keras_api
т
"iter

#beta_1

$beta_2
	%decay
&learning_ratememfmg'mh(mi)mj*mk+ml,mmvnvovp'vq(vr)vs*vt+vu,vv
?
0
'1
(2
)3
*4
+5
,6
7
8
 
?
0
'1
(2
)3
*4
+5
,6
7
8
­
trainable_variables
-metrics

.layers
regularization_losses
	variables
/layer_regularization_losses
0non_trainable_variables
1layer_metrics
 
fd
VARIABLE_VALUEembedding_2/embeddings:layer_with_weights-0/embeddings/.ATTRIBUTES/VARIABLE_VALUE

0
 

0
­
trainable_variables
2metrics

3layers
regularization_losses
	variables
4layer_regularization_losses
5non_trainable_variables
6layer_metrics
~

'kernel
(recurrent_kernel
)bias
7trainable_variables
8regularization_losses
9	variables
:	keras_api
 

'0
(1
)2
 

'0
(1
)2
Й
trainable_variables
;metrics

<states

=layers
regularization_losses
	variables
>layer_regularization_losses
?non_trainable_variables
@layer_metrics
~

*kernel
+recurrent_kernel
,bias
Atrainable_variables
Bregularization_losses
C	variables
D	keras_api
 

*0
+1
,2
 

*0
+1
,2
Й
trainable_variables
Emetrics

Fstates

Glayers
regularization_losses
	variables
Hlayer_regularization_losses
Inon_trainable_variables
Jlayer_metrics
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1
 

0
1
­
trainable_variables
Kmetrics

Llayers
regularization_losses
 	variables
Mlayer_regularization_losses
Nnon_trainable_variables
Olayer_metrics
HF
VARIABLE_VALUE	Adam/iter)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_1+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUE
LJ
VARIABLE_VALUEAdam/beta_2+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUE
JH
VARIABLE_VALUE
Adam/decay*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUE
ZX
VARIABLE_VALUEAdam/learning_rate2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEgru_4/gru_cell_4/kernel0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUE
ge
VARIABLE_VALUE!gru_4/gru_cell_4/recurrent_kernel0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEgru_4/gru_cell_4/bias0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUE
][
VARIABLE_VALUEgru_5/gru_cell_5/kernel0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUE
ge
VARIABLE_VALUE!gru_5/gru_cell_5/recurrent_kernel0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUE
[Y
VARIABLE_VALUEgru_5/gru_cell_5/bias0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUE

P0
Q1

0
1
2
3
 
 
 
 
 
 
 
 

'0
(1
)2
 

'0
(1
)2
­
7trainable_variables
Rmetrics

Slayers
8regularization_losses
9	variables
Tlayer_regularization_losses
Unon_trainable_variables
Vlayer_metrics
 
 

0
 
 
 

*0
+1
,2
 

*0
+1
,2
­
Atrainable_variables
Wmetrics

Xlayers
Bregularization_losses
C	variables
Ylayer_regularization_losses
Znon_trainable_variables
[layer_metrics
 
 

0
 
 
 
 
 
 
 
 
4
	\total
	]count
^	variables
_	keras_api
D
	`total
	acount
b
_fn_kwargs
c	variables
d	keras_api
 
 
 
 
 
 
 
 
 
 
OM
VARIABLE_VALUEtotal4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE
OM
VARIABLE_VALUEcount4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE

\0
]1

^	variables
QO
VARIABLE_VALUEtotal_14keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE
QO
VARIABLE_VALUEcount_14keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE
 

`0
a1

c	variables

VARIABLE_VALUEAdam/embedding_2/embeddings/mVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/mRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/mPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/gru_4/gru_cell_4/kernel/mLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE(Adam/gru_4/gru_cell_4/recurrent_kernel/mLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/gru_4/gru_cell_4/bias/mLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/gru_5/gru_cell_5/kernel/mLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE(Adam/gru_5/gru_cell_5/recurrent_kernel/mLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/gru_5/gru_cell_5/bias/mLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUEAdam/embedding_2/embeddings/vVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
}{
VARIABLE_VALUEAdam/dense_2/kernel/vRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
yw
VARIABLE_VALUEAdam/dense_2/bias/vPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/gru_4/gru_cell_4/kernel/vLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE(Adam/gru_4/gru_cell_4/recurrent_kernel/vLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/gru_4/gru_cell_4/bias/vLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~
VARIABLE_VALUEAdam/gru_5/gru_cell_5/kernel/vLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

VARIABLE_VALUE(Adam/gru_5/gru_cell_5/recurrent_kernel/vLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE
~|
VARIABLE_VALUEAdam/gru_5/gru_cell_5/bias/vLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUE

!serving_default_embedding_2_inputPlaceholder*(
_output_shapes
:џџџџџџџџџш*
dtype0*
shape:џџџџџџџџџш
Љ
StatefulPartitionedCallStatefulPartitionedCall!serving_default_embedding_2_inputembedding_2/embeddingsgru_4/gru_cell_4/biasgru_4/gru_cell_4/kernel!gru_4/gru_cell_4/recurrent_kernelgru_5/gru_cell_5/biasgru_5/gru_cell_5/kernel!gru_5/gru_cell_5/recurrent_kerneldense_2/kerneldense_2/bias*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *,
f'R%
#__inference_signature_wrapper_70250
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 

StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename*embedding_2/embeddings/Read/ReadVariableOp"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOpAdam/iter/Read/ReadVariableOpAdam/beta_1/Read/ReadVariableOpAdam/beta_2/Read/ReadVariableOpAdam/decay/Read/ReadVariableOp&Adam/learning_rate/Read/ReadVariableOp+gru_4/gru_cell_4/kernel/Read/ReadVariableOp5gru_4/gru_cell_4/recurrent_kernel/Read/ReadVariableOp)gru_4/gru_cell_4/bias/Read/ReadVariableOp+gru_5/gru_cell_5/kernel/Read/ReadVariableOp5gru_5/gru_cell_5/recurrent_kernel/Read/ReadVariableOp)gru_5/gru_cell_5/bias/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOp1Adam/embedding_2/embeddings/m/Read/ReadVariableOp)Adam/dense_2/kernel/m/Read/ReadVariableOp'Adam/dense_2/bias/m/Read/ReadVariableOp2Adam/gru_4/gru_cell_4/kernel/m/Read/ReadVariableOp<Adam/gru_4/gru_cell_4/recurrent_kernel/m/Read/ReadVariableOp0Adam/gru_4/gru_cell_4/bias/m/Read/ReadVariableOp2Adam/gru_5/gru_cell_5/kernel/m/Read/ReadVariableOp<Adam/gru_5/gru_cell_5/recurrent_kernel/m/Read/ReadVariableOp0Adam/gru_5/gru_cell_5/bias/m/Read/ReadVariableOp1Adam/embedding_2/embeddings/v/Read/ReadVariableOp)Adam/dense_2/kernel/v/Read/ReadVariableOp'Adam/dense_2/bias/v/Read/ReadVariableOp2Adam/gru_4/gru_cell_4/kernel/v/Read/ReadVariableOp<Adam/gru_4/gru_cell_4/recurrent_kernel/v/Read/ReadVariableOp0Adam/gru_4/gru_cell_4/bias/v/Read/ReadVariableOp2Adam/gru_5/gru_cell_5/kernel/v/Read/ReadVariableOp<Adam/gru_5/gru_cell_5/recurrent_kernel/v/Read/ReadVariableOp0Adam/gru_5/gru_cell_5/bias/v/Read/ReadVariableOpConst*1
Tin*
(2&	*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 *'
f"R 
__inference__traced_save_73142
А	
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameembedding_2/embeddingsdense_2/kerneldense_2/bias	Adam/iterAdam/beta_1Adam/beta_2
Adam/decayAdam/learning_rategru_4/gru_cell_4/kernel!gru_4/gru_cell_4/recurrent_kernelgru_4/gru_cell_4/biasgru_5/gru_cell_5/kernel!gru_5/gru_cell_5/recurrent_kernelgru_5/gru_cell_5/biastotalcounttotal_1count_1Adam/embedding_2/embeddings/mAdam/dense_2/kernel/mAdam/dense_2/bias/mAdam/gru_4/gru_cell_4/kernel/m(Adam/gru_4/gru_cell_4/recurrent_kernel/mAdam/gru_4/gru_cell_4/bias/mAdam/gru_5/gru_cell_5/kernel/m(Adam/gru_5/gru_cell_5/recurrent_kernel/mAdam/gru_5/gru_cell_5/bias/mAdam/embedding_2/embeddings/vAdam/dense_2/kernel/vAdam/dense_2/bias/vAdam/gru_4/gru_cell_4/kernel/v(Adam/gru_4/gru_cell_4/recurrent_kernel/vAdam/gru_4/gru_cell_4/bias/vAdam/gru_5/gru_cell_5/kernel/v(Adam/gru_5/gru_cell_5/recurrent_kernel/vAdam/gru_5/gru_cell_5/bias/v*0
Tin)
'2%*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8 **
f%R#
!__inference__traced_restore_73260Ѕ(
;
щ
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72799

inputs
states_0
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likec
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Const
dropout/MulMulones_like:output:0dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul`
dropout/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout/Shapeг
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2щБо2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/yО
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul_1g
dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_1/Const
dropout_1/MulMulones_like:output:0dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Muld
dropout_1/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_1/Shapeи
&dropout_1/random_uniform/RandomUniformRandomUniformdropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2}2(
&dropout_1/random_uniform/RandomUniformy
dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_1/GreaterEqual/yЦ
dropout_1/GreaterEqualGreaterEqual/dropout_1/random_uniform/RandomUniform:output:0!dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/GreaterEqual
dropout_1/CastCastdropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Cast
dropout_1/Mul_1Muldropout_1/Mul:z:0dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Mul_1g
dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_2/Const
dropout_2/MulMulones_like:output:0dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Muld
dropout_2/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_2/Shapeй
&dropout_2/random_uniform/RandomUniformRandomUniformdropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Їќ2(
&dropout_2/random_uniform/RandomUniformy
dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_2/GreaterEqual/yЦ
dropout_2/GreaterEqualGreaterEqual/dropout_2/random_uniform/RandomUniform:output:0!dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/GreaterEqual
dropout_2/CastCastdropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Cast
dropout_2/Mul_1Muldropout_2/Mul:z:0dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Mul_1x
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack^
mulMulinputsdropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOp{
MatMul_1MatMulstates_0MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2^
mul_2MulSigmoid:y:0states_0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
Нк
З
G__inference_sequential_2_layer_call_and_return_conditional_losses_71032

inputs&
"embedding_2_embedding_lookup_70693,
(gru_4_gru_cell_4_readvariableop_resource3
/gru_4_gru_cell_4_matmul_readvariableop_resource5
1gru_4_gru_cell_4_matmul_1_readvariableop_resource,
(gru_5_gru_cell_5_readvariableop_resource3
/gru_5_gru_cell_5_matmul_readvariableop_resource5
1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource
identityЂgru_4/whileЂgru_5/whilev
embedding_2/CastCastinputs*

DstT0*

SrcT0*(
_output_shapes
:џџџџџџџџџш2
embedding_2/Cast
embedding_2/embedding_lookupResourceGather"embedding_2_embedding_lookup_70693embedding_2/Cast:y:0*
Tindices0*5
_class+
)'loc:@embedding_2/embedding_lookup/70693*,
_output_shapes
:џџџџџџџџџш*
dtype02
embedding_2/embedding_lookupя
%embedding_2/embedding_lookup/IdentityIdentity%embedding_2/embedding_lookup:output:0*
T0*5
_class+
)'loc:@embedding_2/embedding_lookup/70693*,
_output_shapes
:џџџџџџџџџш2'
%embedding_2/embedding_lookup/IdentityХ
'embedding_2/embedding_lookup/Identity_1Identity.embedding_2/embedding_lookup/Identity:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2)
'embedding_2/embedding_lookup/Identity_1z
gru_4/ShapeShape0embedding_2/embedding_lookup/Identity_1:output:0*
T0*
_output_shapes
:2
gru_4/Shape
gru_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice/stack
gru_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice/stack_1
gru_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice/stack_2
gru_4/strided_sliceStridedSlicegru_4/Shape:output:0"gru_4/strided_slice/stack:output:0$gru_4/strided_slice/stack_1:output:0$gru_4/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_4/strided_sliceh
gru_4/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/zeros/mul/y
gru_4/zeros/mulMulgru_4/strided_slice:output:0gru_4/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
gru_4/zeros/mulk
gru_4/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
gru_4/zeros/Less/y
gru_4/zeros/LessLessgru_4/zeros/mul:z:0gru_4/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
gru_4/zeros/Lessn
gru_4/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
gru_4/zeros/packed/1
gru_4/zeros/packedPackgru_4/strided_slice:output:0gru_4/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
gru_4/zeros/packedk
gru_4/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_4/zeros/Const
gru_4/zerosFillgru_4/zeros/packed:output:0gru_4/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/zeros
gru_4/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_4/transpose/permЗ
gru_4/transpose	Transpose0embedding_2/embedding_lookup/Identity_1:output:0gru_4/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
gru_4/transposea
gru_4/Shape_1Shapegru_4/transpose:y:0*
T0*
_output_shapes
:2
gru_4/Shape_1
gru_4/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_1/stack
gru_4/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_1/stack_1
gru_4/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_1/stack_2
gru_4/strided_slice_1StridedSlicegru_4/Shape_1:output:0$gru_4/strided_slice_1/stack:output:0&gru_4/strided_slice_1/stack_1:output:0&gru_4/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_4/strided_slice_1
!gru_4/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2#
!gru_4/TensorArrayV2/element_shapeЪ
gru_4/TensorArrayV2TensorListReserve*gru_4/TensorArrayV2/element_shape:output:0gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_4/TensorArrayV2Ы
;gru_4/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2=
;gru_4/TensorArrayUnstack/TensorListFromTensor/element_shape
-gru_4/TensorArrayUnstack/TensorListFromTensorTensorListFromTensorgru_4/transpose:y:0Dgru_4/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02/
-gru_4/TensorArrayUnstack/TensorListFromTensor
gru_4/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_2/stack
gru_4/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_2/stack_1
gru_4/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_2/stack_2 
gru_4/strided_slice_2StridedSlicegru_4/transpose:y:0$gru_4/strided_slice_2/stack:output:0&gru_4/strided_slice_2/stack_1:output:0&gru_4/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_4/strided_slice_2
 gru_4/gru_cell_4/ones_like/ShapeShapegru_4/strided_slice_2:output:0*
T0*
_output_shapes
:2"
 gru_4/gru_cell_4/ones_like/Shape
 gru_4/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 gru_4/gru_cell_4/ones_like/ConstШ
gru_4/gru_cell_4/ones_likeFill)gru_4/gru_cell_4/ones_like/Shape:output:0)gru_4/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/ones_likeЋ
gru_4/gru_cell_4/ReadVariableOpReadVariableOp(gru_4_gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02!
gru_4/gru_cell_4/ReadVariableOp
gru_4/gru_cell_4/unstackUnpack'gru_4/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_4/gru_cell_4/unstackЊ
gru_4/gru_cell_4/mulMulgru_4/strided_slice_2:output:0#gru_4/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mulР
&gru_4/gru_cell_4/MatMul/ReadVariableOpReadVariableOp/gru_4_gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02(
&gru_4/gru_cell_4/MatMul/ReadVariableOpИ
gru_4/gru_cell_4/MatMulMatMulgru_4/gru_cell_4/mul:z:0.gru_4/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/MatMulЗ
gru_4/gru_cell_4/BiasAddBiasAdd!gru_4/gru_cell_4/MatMul:product:0!gru_4/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/BiasAddr
gru_4/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/gru_cell_4/Const
 gru_4/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 gru_4/gru_cell_4/split/split_dim№
gru_4/gru_cell_4/splitSplit)gru_4/gru_cell_4/split/split_dim:output:0!gru_4/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/gru_cell_4/splitЦ
(gru_4/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp1gru_4_gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02*
(gru_4/gru_cell_4/MatMul_1/ReadVariableOpК
gru_4/gru_cell_4/MatMul_1MatMulgru_4/zeros:output:00gru_4/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/MatMul_1Н
gru_4/gru_cell_4/BiasAdd_1BiasAdd#gru_4/gru_cell_4/MatMul_1:product:0!gru_4/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/BiasAdd_1
gru_4/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_4/gru_cell_4/Const_1
"gru_4/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"gru_4/gru_cell_4/split_1/split_dimЈ
gru_4/gru_cell_4/split_1SplitV#gru_4/gru_cell_4/BiasAdd_1:output:0!gru_4/gru_cell_4/Const_1:output:0+gru_4/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/gru_cell_4/split_1Ћ
gru_4/gru_cell_4/addAddV2gru_4/gru_cell_4/split:output:0!gru_4/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add
gru_4/gru_cell_4/SigmoidSigmoidgru_4/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/SigmoidЏ
gru_4/gru_cell_4/add_1AddV2gru_4/gru_cell_4/split:output:1!gru_4/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_1
gru_4/gru_cell_4/Sigmoid_1Sigmoidgru_4/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/Sigmoid_1Ќ
gru_4/gru_cell_4/mul_1Mulgru_4/gru_cell_4/Sigmoid_1:y:0!gru_4/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_1Ј
gru_4/gru_cell_4/add_2AddV2gru_4/gru_cell_4/split:output:2gru_4/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_2
gru_4/gru_cell_4/Sigmoid_2Sigmoidgru_4/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/Sigmoid_2
gru_4/gru_cell_4/mul_2Mulgru_4/gru_cell_4/Sigmoid:y:0gru_4/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_2u
gru_4/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_4/gru_cell_4/sub/xЄ
gru_4/gru_cell_4/subSubgru_4/gru_cell_4/sub/x:output:0gru_4/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/subЃ
gru_4/gru_cell_4/mul_3Mulgru_4/gru_cell_4/sub:z:0gru_4/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_3Ѓ
gru_4/gru_cell_4/add_3AddV2gru_4/gru_cell_4/mul_2:z:0gru_4/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_3
#gru_4/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2%
#gru_4/TensorArrayV2_1/element_shapeа
gru_4/TensorArrayV2_1TensorListReserve,gru_4/TensorArrayV2_1/element_shape:output:0gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_4/TensorArrayV2_1Z

gru_4/timeConst*
_output_shapes
: *
dtype0*
value	B : 2

gru_4/time
gru_4/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2 
gru_4/while/maximum_iterationsv
gru_4/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
gru_4/while/loop_counterљ
gru_4/whileWhile!gru_4/while/loop_counter:output:0'gru_4/while/maximum_iterations:output:0gru_4/time:output:0gru_4/TensorArrayV2_1:handle:0gru_4/zeros:output:0gru_4/strided_slice_1:output:0=gru_4/TensorArrayUnstack/TensorListFromTensor:output_handle:0(gru_4_gru_cell_4_readvariableop_resource/gru_4_gru_cell_4_matmul_readvariableop_resource1gru_4_gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*"
bodyR
gru_4_while_body_70768*"
condR
gru_4_while_cond_70767*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
gru_4/whileС
6gru_4/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   28
6gru_4/TensorArrayV2Stack/TensorListStack/element_shape
(gru_4/TensorArrayV2Stack/TensorListStackTensorListStackgru_4/while:output:3?gru_4/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02*
(gru_4/TensorArrayV2Stack/TensorListStack
gru_4/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
gru_4/strided_slice_3/stack
gru_4/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_3/stack_1
gru_4/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_3/stack_2О
gru_4/strided_slice_3StridedSlice1gru_4/TensorArrayV2Stack/TensorListStack:tensor:0$gru_4/strided_slice_3/stack:output:0&gru_4/strided_slice_3/stack_1:output:0&gru_4/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_4/strided_slice_3
gru_4/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_4/transpose_1/permО
gru_4/transpose_1	Transpose1gru_4/TensorArrayV2Stack/TensorListStack:tensor:0gru_4/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
gru_4/transpose_1r
gru_4/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_4/runtime_
gru_5/ShapeShapegru_4/transpose_1:y:0*
T0*
_output_shapes
:2
gru_5/Shape
gru_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice/stack
gru_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice/stack_1
gru_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice/stack_2
gru_5/strided_sliceStridedSlicegru_5/Shape:output:0"gru_5/strided_slice/stack:output:0$gru_5/strided_slice/stack_1:output:0$gru_5/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_5/strided_sliceh
gru_5/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/zeros/mul/y
gru_5/zeros/mulMulgru_5/strided_slice:output:0gru_5/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
gru_5/zeros/mulk
gru_5/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
gru_5/zeros/Less/y
gru_5/zeros/LessLessgru_5/zeros/mul:z:0gru_5/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
gru_5/zeros/Lessn
gru_5/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
gru_5/zeros/packed/1
gru_5/zeros/packedPackgru_5/strided_slice:output:0gru_5/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
gru_5/zeros/packedk
gru_5/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_5/zeros/Const
gru_5/zerosFillgru_5/zeros/packed:output:0gru_5/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/zeros
gru_5/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_5/transpose/perm
gru_5/transpose	Transposegru_4/transpose_1:y:0gru_5/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
gru_5/transposea
gru_5/Shape_1Shapegru_5/transpose:y:0*
T0*
_output_shapes
:2
gru_5/Shape_1
gru_5/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_1/stack
gru_5/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_1/stack_1
gru_5/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_1/stack_2
gru_5/strided_slice_1StridedSlicegru_5/Shape_1:output:0$gru_5/strided_slice_1/stack:output:0&gru_5/strided_slice_1/stack_1:output:0&gru_5/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_5/strided_slice_1
!gru_5/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2#
!gru_5/TensorArrayV2/element_shapeЪ
gru_5/TensorArrayV2TensorListReserve*gru_5/TensorArrayV2/element_shape:output:0gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_5/TensorArrayV2Ы
;gru_5/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2=
;gru_5/TensorArrayUnstack/TensorListFromTensor/element_shape
-gru_5/TensorArrayUnstack/TensorListFromTensorTensorListFromTensorgru_5/transpose:y:0Dgru_5/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02/
-gru_5/TensorArrayUnstack/TensorListFromTensor
gru_5/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_2/stack
gru_5/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_2/stack_1
gru_5/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_2/stack_2 
gru_5/strided_slice_2StridedSlicegru_5/transpose:y:0$gru_5/strided_slice_2/stack:output:0&gru_5/strided_slice_2/stack_1:output:0&gru_5/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_5/strided_slice_2
 gru_5/gru_cell_5/ones_like/ShapeShapegru_5/strided_slice_2:output:0*
T0*
_output_shapes
:2"
 gru_5/gru_cell_5/ones_like/Shape
 gru_5/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 gru_5/gru_cell_5/ones_like/ConstШ
gru_5/gru_cell_5/ones_likeFill)gru_5/gru_cell_5/ones_like/Shape:output:0)gru_5/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/ones_likeЋ
gru_5/gru_cell_5/ReadVariableOpReadVariableOp(gru_5_gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02!
gru_5/gru_cell_5/ReadVariableOp
gru_5/gru_cell_5/unstackUnpack'gru_5/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_5/gru_cell_5/unstackЊ
gru_5/gru_cell_5/mulMulgru_5/strided_slice_2:output:0#gru_5/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mulР
&gru_5/gru_cell_5/MatMul/ReadVariableOpReadVariableOp/gru_5_gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02(
&gru_5/gru_cell_5/MatMul/ReadVariableOpИ
gru_5/gru_cell_5/MatMulMatMulgru_5/gru_cell_5/mul:z:0.gru_5/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/MatMulЗ
gru_5/gru_cell_5/BiasAddBiasAdd!gru_5/gru_cell_5/MatMul:product:0!gru_5/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/BiasAddr
gru_5/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/gru_cell_5/Const
 gru_5/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 gru_5/gru_cell_5/split/split_dim№
gru_5/gru_cell_5/splitSplit)gru_5/gru_cell_5/split/split_dim:output:0!gru_5/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/gru_cell_5/splitЦ
(gru_5/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02*
(gru_5/gru_cell_5/MatMul_1/ReadVariableOpК
gru_5/gru_cell_5/MatMul_1MatMulgru_5/zeros:output:00gru_5/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/MatMul_1Н
gru_5/gru_cell_5/BiasAdd_1BiasAdd#gru_5/gru_cell_5/MatMul_1:product:0!gru_5/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/BiasAdd_1
gru_5/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_5/gru_cell_5/Const_1
"gru_5/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"gru_5/gru_cell_5/split_1/split_dimЈ
gru_5/gru_cell_5/split_1SplitV#gru_5/gru_cell_5/BiasAdd_1:output:0!gru_5/gru_cell_5/Const_1:output:0+gru_5/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/gru_cell_5/split_1Ћ
gru_5/gru_cell_5/addAddV2gru_5/gru_cell_5/split:output:0!gru_5/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add
gru_5/gru_cell_5/SigmoidSigmoidgru_5/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/SigmoidЏ
gru_5/gru_cell_5/add_1AddV2gru_5/gru_cell_5/split:output:1!gru_5/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_1
gru_5/gru_cell_5/Sigmoid_1Sigmoidgru_5/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/Sigmoid_1Ќ
gru_5/gru_cell_5/mul_1Mulgru_5/gru_cell_5/Sigmoid_1:y:0!gru_5/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_1Ј
gru_5/gru_cell_5/add_2AddV2gru_5/gru_cell_5/split:output:2gru_5/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_2
gru_5/gru_cell_5/Sigmoid_2Sigmoidgru_5/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/Sigmoid_2
gru_5/gru_cell_5/mul_2Mulgru_5/gru_cell_5/Sigmoid:y:0gru_5/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_2u
gru_5/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_5/gru_cell_5/sub/xЄ
gru_5/gru_cell_5/subSubgru_5/gru_cell_5/sub/x:output:0gru_5/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/subЃ
gru_5/gru_cell_5/mul_3Mulgru_5/gru_cell_5/sub:z:0gru_5/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_3Ѓ
gru_5/gru_cell_5/add_3AddV2gru_5/gru_cell_5/mul_2:z:0gru_5/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_3
#gru_5/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2%
#gru_5/TensorArrayV2_1/element_shapeа
gru_5/TensorArrayV2_1TensorListReserve,gru_5/TensorArrayV2_1/element_shape:output:0gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_5/TensorArrayV2_1Z

gru_5/timeConst*
_output_shapes
: *
dtype0*
value	B : 2

gru_5/time
gru_5/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2 
gru_5/while/maximum_iterationsv
gru_5/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
gru_5/while/loop_counterљ
gru_5/whileWhile!gru_5/while/loop_counter:output:0'gru_5/while/maximum_iterations:output:0gru_5/time:output:0gru_5/TensorArrayV2_1:handle:0gru_5/zeros:output:0gru_5/strided_slice_1:output:0=gru_5/TensorArrayUnstack/TensorListFromTensor:output_handle:0(gru_5_gru_cell_5_readvariableop_resource/gru_5_gru_cell_5_matmul_readvariableop_resource1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*"
bodyR
gru_5_while_body_70931*"
condR
gru_5_while_cond_70930*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
gru_5/whileС
6gru_5/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   28
6gru_5/TensorArrayV2Stack/TensorListStack/element_shape
(gru_5/TensorArrayV2Stack/TensorListStackTensorListStackgru_5/while:output:3?gru_5/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02*
(gru_5/TensorArrayV2Stack/TensorListStack
gru_5/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
gru_5/strided_slice_3/stack
gru_5/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_3/stack_1
gru_5/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_3/stack_2О
gru_5/strided_slice_3StridedSlice1gru_5/TensorArrayV2Stack/TensorListStack:tensor:0$gru_5/strided_slice_3/stack:output:0&gru_5/strided_slice_3/stack_1:output:0&gru_5/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_5/strided_slice_3
gru_5/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_5/transpose_1/permО
gru_5/transpose_1	Transpose1gru_5/TensorArrayV2Stack/TensorListStack:tensor:0gru_5/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
gru_5/transpose_1r
gru_5/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_5/runtimeЅ
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_2/MatMul/ReadVariableOpЃ
dense_2/MatMulMatMulgru_5/strided_slice_3:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/MatMulЄ
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_2/BiasAdd/ReadVariableOpЁ
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/BiasAddy
dense_2/SigmoidSigmoiddense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/Sigmoid
IdentityIdentitydense_2/Sigmoid:y:0^gru_4/while^gru_5/while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2
gru_4/whilegru_4/while2
gru_5/whilegru_5/while:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ж
ь
#__inference_signature_wrapper_70250
embedding_2_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
identityЂStatefulPartitionedCallЖ
StatefulPartitionedCallStatefulPartitionedCallembedding_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *)
f$R"
 __inference__wrapped_model_680282
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input


%__inference_gru_4_layer_call_fn_71488

inputs
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_694562
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
;
ч
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_68128

inputs

states
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likec
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Const
dropout/MulMulones_like:output:0dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul`
dropout/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout/Shapeв
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2туB2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/yО
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul_1g
dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_1/Const
dropout_1/MulMulones_like:output:0dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Muld
dropout_1/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_1/Shapeй
&dropout_1/random_uniform/RandomUniformRandomUniformdropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2хщМ2(
&dropout_1/random_uniform/RandomUniformy
dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_1/GreaterEqual/yЦ
dropout_1/GreaterEqualGreaterEqual/dropout_1/random_uniform/RandomUniform:output:0!dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/GreaterEqual
dropout_1/CastCastdropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Cast
dropout_1/Mul_1Muldropout_1/Mul:z:0dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Mul_1g
dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_2/Const
dropout_2/MulMulones_like:output:0dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Muld
dropout_2/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_2/Shapeй
&dropout_2/random_uniform/RandomUniformRandomUniformdropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ыцм2(
&dropout_2/random_uniform/RandomUniformy
dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_2/GreaterEqual/yЦ
dropout_2/GreaterEqualGreaterEqual/dropout_2/random_uniform/RandomUniform:output:0!dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/GreaterEqual
dropout_2/CastCastdropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Cast
dropout_2/Mul_1Muldropout_2/Mul:z:0dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Mul_1x
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack^
mulMulinputsdropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOpy
MatMul_1MatMulstatesMatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2\
mul_2MulSigmoid:y:0states*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:OK
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_namestates
гi
Ў
while_body_72404
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like
while/gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_5/dropout/ConstУ
while/gru_cell_5/dropout/MulMul#while/gru_cell_5/ones_like:output:0'while/gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/Mul
while/gru_cell_5/dropout/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_5/dropout/Shape
5while/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2пъ27
5while/gru_cell_5/dropout/random_uniform/RandomUniform
'while/gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_5/dropout/GreaterEqual/y
%while/gru_cell_5/dropout/GreaterEqualGreaterEqual>while/gru_cell_5/dropout/random_uniform/RandomUniform:output:00while/gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_5/dropout/GreaterEqualВ
while/gru_cell_5/dropout/CastCast)while/gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/CastО
while/gru_cell_5/dropout/Mul_1Mul while/gru_cell_5/dropout/Mul:z:0!while/gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout/Mul_1
 while/gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_1/ConstЩ
while/gru_cell_5/dropout_1/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_1/Mul
 while/gru_cell_5/dropout_1/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_1/Shape
7while/gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Ш29
7while/gru_cell_5/dropout_1/random_uniform/RandomUniform
)while/gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_1/GreaterEqual/y
'while/gru_cell_5/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_1/GreaterEqualИ
while/gru_cell_5/dropout_1/CastCast+while/gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_1/CastЦ
 while/gru_cell_5/dropout_1/Mul_1Mul"while/gru_cell_5/dropout_1/Mul:z:0#while/gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_1/Mul_1
 while/gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_2/ConstЩ
while/gru_cell_5/dropout_2/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_2/Mul
 while/gru_cell_5/dropout_2/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_2/Shape
7while/gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Рћ29
7while/gru_cell_5/dropout_2/random_uniform/RandomUniform
)while/gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_2/GreaterEqual/y
'while/gru_cell_5/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_2/GreaterEqualИ
while/gru_cell_5/dropout_2/CastCast+while/gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_2/CastЦ
 while/gru_cell_5/dropout_2/Mul_1Mul"while/gru_cell_5/dropout_2/Mul:z:0#while/gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_2/Mul_1­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackЛ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
о
Њ
 __inference__wrapped_model_68028
embedding_2_input3
/sequential_2_embedding_2_embedding_lookup_676899
5sequential_2_gru_4_gru_cell_4_readvariableop_resource@
<sequential_2_gru_4_gru_cell_4_matmul_readvariableop_resourceB
>sequential_2_gru_4_gru_cell_4_matmul_1_readvariableop_resource9
5sequential_2_gru_5_gru_cell_5_readvariableop_resource@
<sequential_2_gru_5_gru_cell_5_matmul_readvariableop_resourceB
>sequential_2_gru_5_gru_cell_5_matmul_1_readvariableop_resource7
3sequential_2_dense_2_matmul_readvariableop_resource8
4sequential_2_dense_2_biasadd_readvariableop_resource
identityЂsequential_2/gru_4/whileЂsequential_2/gru_5/while
sequential_2/embedding_2/CastCastembedding_2_input*

DstT0*

SrcT0*(
_output_shapes
:џџџџџџџџџш2
sequential_2/embedding_2/CastЫ
)sequential_2/embedding_2/embedding_lookupResourceGather/sequential_2_embedding_2_embedding_lookup_67689!sequential_2/embedding_2/Cast:y:0*
Tindices0*B
_class8
64loc:@sequential_2/embedding_2/embedding_lookup/67689*,
_output_shapes
:џџџџџџџџџш*
dtype02+
)sequential_2/embedding_2/embedding_lookupЃ
2sequential_2/embedding_2/embedding_lookup/IdentityIdentity2sequential_2/embedding_2/embedding_lookup:output:0*
T0*B
_class8
64loc:@sequential_2/embedding_2/embedding_lookup/67689*,
_output_shapes
:џџџџџџџџџш24
2sequential_2/embedding_2/embedding_lookup/Identityь
4sequential_2/embedding_2/embedding_lookup/Identity_1Identity;sequential_2/embedding_2/embedding_lookup/Identity:output:0*
T0*,
_output_shapes
:џџџџџџџџџш26
4sequential_2/embedding_2/embedding_lookup/Identity_1Ё
sequential_2/gru_4/ShapeShape=sequential_2/embedding_2/embedding_lookup/Identity_1:output:0*
T0*
_output_shapes
:2
sequential_2/gru_4/Shape
&sequential_2/gru_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&sequential_2/gru_4/strided_slice/stack
(sequential_2/gru_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(sequential_2/gru_4/strided_slice/stack_1
(sequential_2/gru_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(sequential_2/gru_4/strided_slice/stack_2д
 sequential_2/gru_4/strided_sliceStridedSlice!sequential_2/gru_4/Shape:output:0/sequential_2/gru_4/strided_slice/stack:output:01sequential_2/gru_4/strided_slice/stack_1:output:01sequential_2/gru_4/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 sequential_2/gru_4/strided_slice
sequential_2/gru_4/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2 
sequential_2/gru_4/zeros/mul/yИ
sequential_2/gru_4/zeros/mulMul)sequential_2/gru_4/strided_slice:output:0'sequential_2/gru_4/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_4/zeros/mul
sequential_2/gru_4/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2!
sequential_2/gru_4/zeros/Less/yГ
sequential_2/gru_4/zeros/LessLess sequential_2/gru_4/zeros/mul:z:0(sequential_2/gru_4/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_4/zeros/Less
!sequential_2/gru_4/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2#
!sequential_2/gru_4/zeros/packed/1Я
sequential_2/gru_4/zeros/packedPack)sequential_2/gru_4/strided_slice:output:0*sequential_2/gru_4/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2!
sequential_2/gru_4/zeros/packed
sequential_2/gru_4/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
sequential_2/gru_4/zeros/ConstС
sequential_2/gru_4/zerosFill(sequential_2/gru_4/zeros/packed:output:0'sequential_2/gru_4/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sequential_2/gru_4/zeros
!sequential_2/gru_4/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2#
!sequential_2/gru_4/transpose/permы
sequential_2/gru_4/transpose	Transpose=sequential_2/embedding_2/embedding_lookup/Identity_1:output:0*sequential_2/gru_4/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
sequential_2/gru_4/transpose
sequential_2/gru_4/Shape_1Shape sequential_2/gru_4/transpose:y:0*
T0*
_output_shapes
:2
sequential_2/gru_4/Shape_1
(sequential_2/gru_4/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(sequential_2/gru_4/strided_slice_1/stackЂ
*sequential_2/gru_4/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_4/strided_slice_1/stack_1Ђ
*sequential_2/gru_4/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_4/strided_slice_1/stack_2р
"sequential_2/gru_4/strided_slice_1StridedSlice#sequential_2/gru_4/Shape_1:output:01sequential_2/gru_4/strided_slice_1/stack:output:03sequential_2/gru_4/strided_slice_1/stack_1:output:03sequential_2/gru_4/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"sequential_2/gru_4/strided_slice_1Ћ
.sequential_2/gru_4/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ20
.sequential_2/gru_4/TensorArrayV2/element_shapeў
 sequential_2/gru_4/TensorArrayV2TensorListReserve7sequential_2/gru_4/TensorArrayV2/element_shape:output:0+sequential_2/gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02"
 sequential_2/gru_4/TensorArrayV2х
Hsequential_2/gru_4/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2J
Hsequential_2/gru_4/TensorArrayUnstack/TensorListFromTensor/element_shapeФ
:sequential_2/gru_4/TensorArrayUnstack/TensorListFromTensorTensorListFromTensor sequential_2/gru_4/transpose:y:0Qsequential_2/gru_4/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02<
:sequential_2/gru_4/TensorArrayUnstack/TensorListFromTensor
(sequential_2/gru_4/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(sequential_2/gru_4/strided_slice_2/stackЂ
*sequential_2/gru_4/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_4/strided_slice_2/stack_1Ђ
*sequential_2/gru_4/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_4/strided_slice_2/stack_2ю
"sequential_2/gru_4/strided_slice_2StridedSlice sequential_2/gru_4/transpose:y:01sequential_2/gru_4/strided_slice_2/stack:output:03sequential_2/gru_4/strided_slice_2/stack_1:output:03sequential_2/gru_4/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2$
"sequential_2/gru_4/strided_slice_2Й
-sequential_2/gru_4/gru_cell_4/ones_like/ShapeShape+sequential_2/gru_4/strided_slice_2:output:0*
T0*
_output_shapes
:2/
-sequential_2/gru_4/gru_cell_4/ones_like/ShapeЃ
-sequential_2/gru_4/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2/
-sequential_2/gru_4/gru_cell_4/ones_like/Constќ
'sequential_2/gru_4/gru_cell_4/ones_likeFill6sequential_2/gru_4/gru_cell_4/ones_like/Shape:output:06sequential_2/gru_4/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/gru_cell_4/ones_likeв
,sequential_2/gru_4/gru_cell_4/ReadVariableOpReadVariableOp5sequential_2_gru_4_gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02.
,sequential_2/gru_4/gru_cell_4/ReadVariableOpФ
%sequential_2/gru_4/gru_cell_4/unstackUnpack4sequential_2/gru_4/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2'
%sequential_2/gru_4/gru_cell_4/unstackо
!sequential_2/gru_4/gru_cell_4/mulMul+sequential_2/gru_4/strided_slice_2:output:00sequential_2/gru_4/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_4/gru_cell_4/mulч
3sequential_2/gru_4/gru_cell_4/MatMul/ReadVariableOpReadVariableOp<sequential_2_gru_4_gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype025
3sequential_2/gru_4/gru_cell_4/MatMul/ReadVariableOpь
$sequential_2/gru_4/gru_cell_4/MatMulMatMul%sequential_2/gru_4/gru_cell_4/mul:z:0;sequential_2/gru_4/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02&
$sequential_2/gru_4/gru_cell_4/MatMulы
%sequential_2/gru_4/gru_cell_4/BiasAddBiasAdd.sequential_2/gru_4/gru_cell_4/MatMul:product:0.sequential_2/gru_4/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02'
%sequential_2/gru_4/gru_cell_4/BiasAdd
#sequential_2/gru_4/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2%
#sequential_2/gru_4/gru_cell_4/ConstЉ
-sequential_2/gru_4/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2/
-sequential_2/gru_4/gru_cell_4/split/split_dimЄ
#sequential_2/gru_4/gru_cell_4/splitSplit6sequential_2/gru_4/gru_cell_4/split/split_dim:output:0.sequential_2/gru_4/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2%
#sequential_2/gru_4/gru_cell_4/splitэ
5sequential_2/gru_4/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp>sequential_2_gru_4_gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype027
5sequential_2/gru_4/gru_cell_4/MatMul_1/ReadVariableOpю
&sequential_2/gru_4/gru_cell_4/MatMul_1MatMul!sequential_2/gru_4/zeros:output:0=sequential_2/gru_4/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02(
&sequential_2/gru_4/gru_cell_4/MatMul_1ё
'sequential_2/gru_4/gru_cell_4/BiasAdd_1BiasAdd0sequential_2/gru_4/gru_cell_4/MatMul_1:product:0.sequential_2/gru_4/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02)
'sequential_2/gru_4/gru_cell_4/BiasAdd_1Ѓ
%sequential_2/gru_4/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2'
%sequential_2/gru_4/gru_cell_4/Const_1­
/sequential_2/gru_4/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ21
/sequential_2/gru_4/gru_cell_4/split_1/split_dimщ
%sequential_2/gru_4/gru_cell_4/split_1SplitV0sequential_2/gru_4/gru_cell_4/BiasAdd_1:output:0.sequential_2/gru_4/gru_cell_4/Const_1:output:08sequential_2/gru_4/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2'
%sequential_2/gru_4/gru_cell_4/split_1п
!sequential_2/gru_4/gru_cell_4/addAddV2,sequential_2/gru_4/gru_cell_4/split:output:0.sequential_2/gru_4/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_4/gru_cell_4/addВ
%sequential_2/gru_4/gru_cell_4/SigmoidSigmoid%sequential_2/gru_4/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%sequential_2/gru_4/gru_cell_4/Sigmoidу
#sequential_2/gru_4/gru_cell_4/add_1AddV2,sequential_2/gru_4/gru_cell_4/split:output:1.sequential_2/gru_4/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/add_1И
'sequential_2/gru_4/gru_cell_4/Sigmoid_1Sigmoid'sequential_2/gru_4/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/gru_cell_4/Sigmoid_1р
#sequential_2/gru_4/gru_cell_4/mul_1Mul+sequential_2/gru_4/gru_cell_4/Sigmoid_1:y:0.sequential_2/gru_4/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/mul_1м
#sequential_2/gru_4/gru_cell_4/add_2AddV2,sequential_2/gru_4/gru_cell_4/split:output:2'sequential_2/gru_4/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/add_2И
'sequential_2/gru_4/gru_cell_4/Sigmoid_2Sigmoid'sequential_2/gru_4/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/gru_cell_4/Sigmoid_2б
#sequential_2/gru_4/gru_cell_4/mul_2Mul)sequential_2/gru_4/gru_cell_4/Sigmoid:y:0!sequential_2/gru_4/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/mul_2
#sequential_2/gru_4/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2%
#sequential_2/gru_4/gru_cell_4/sub/xи
!sequential_2/gru_4/gru_cell_4/subSub,sequential_2/gru_4/gru_cell_4/sub/x:output:0)sequential_2/gru_4/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_4/gru_cell_4/subз
#sequential_2/gru_4/gru_cell_4/mul_3Mul%sequential_2/gru_4/gru_cell_4/sub:z:0+sequential_2/gru_4/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/mul_3з
#sequential_2/gru_4/gru_cell_4/add_3AddV2'sequential_2/gru_4/gru_cell_4/mul_2:z:0'sequential_2/gru_4/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/gru_cell_4/add_3Е
0sequential_2/gru_4/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0sequential_2/gru_4/TensorArrayV2_1/element_shape
"sequential_2/gru_4/TensorArrayV2_1TensorListReserve9sequential_2/gru_4/TensorArrayV2_1/element_shape:output:0+sequential_2/gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02$
"sequential_2/gru_4/TensorArrayV2_1t
sequential_2/gru_4/timeConst*
_output_shapes
: *
dtype0*
value	B : 2
sequential_2/gru_4/timeЅ
+sequential_2/gru_4/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2-
+sequential_2/gru_4/while/maximum_iterations
%sequential_2/gru_4/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2'
%sequential_2/gru_4/while/loop_counterЏ
sequential_2/gru_4/whileWhile.sequential_2/gru_4/while/loop_counter:output:04sequential_2/gru_4/while/maximum_iterations:output:0 sequential_2/gru_4/time:output:0+sequential_2/gru_4/TensorArrayV2_1:handle:0!sequential_2/gru_4/zeros:output:0+sequential_2/gru_4/strided_slice_1:output:0Jsequential_2/gru_4/TensorArrayUnstack/TensorListFromTensor:output_handle:05sequential_2_gru_4_gru_cell_4_readvariableop_resource<sequential_2_gru_4_gru_cell_4_matmul_readvariableop_resource>sequential_2_gru_4_gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*/
body'R%
#sequential_2_gru_4_while_body_67764*/
cond'R%
#sequential_2_gru_4_while_cond_67763*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
sequential_2/gru_4/whileл
Csequential_2/gru_4/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2E
Csequential_2/gru_4/TensorArrayV2Stack/TensorListStack/element_shapeЕ
5sequential_2/gru_4/TensorArrayV2Stack/TensorListStackTensorListStack!sequential_2/gru_4/while:output:3Lsequential_2/gru_4/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype027
5sequential_2/gru_4/TensorArrayV2Stack/TensorListStackЇ
(sequential_2/gru_4/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2*
(sequential_2/gru_4/strided_slice_3/stackЂ
*sequential_2/gru_4/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2,
*sequential_2/gru_4/strided_slice_3/stack_1Ђ
*sequential_2/gru_4/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_4/strided_slice_3/stack_2
"sequential_2/gru_4/strided_slice_3StridedSlice>sequential_2/gru_4/TensorArrayV2Stack/TensorListStack:tensor:01sequential_2/gru_4/strided_slice_3/stack:output:03sequential_2/gru_4/strided_slice_3/stack_1:output:03sequential_2/gru_4/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2$
"sequential_2/gru_4/strided_slice_3
#sequential_2/gru_4/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2%
#sequential_2/gru_4/transpose_1/permђ
sequential_2/gru_4/transpose_1	Transpose>sequential_2/gru_4/TensorArrayV2Stack/TensorListStack:tensor:0,sequential_2/gru_4/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2 
sequential_2/gru_4/transpose_1
sequential_2/gru_4/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
sequential_2/gru_4/runtime
sequential_2/gru_5/ShapeShape"sequential_2/gru_4/transpose_1:y:0*
T0*
_output_shapes
:2
sequential_2/gru_5/Shape
&sequential_2/gru_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2(
&sequential_2/gru_5/strided_slice/stack
(sequential_2/gru_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2*
(sequential_2/gru_5/strided_slice/stack_1
(sequential_2/gru_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2*
(sequential_2/gru_5/strided_slice/stack_2д
 sequential_2/gru_5/strided_sliceStridedSlice!sequential_2/gru_5/Shape:output:0/sequential_2/gru_5/strided_slice/stack:output:01sequential_2/gru_5/strided_slice/stack_1:output:01sequential_2/gru_5/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2"
 sequential_2/gru_5/strided_slice
sequential_2/gru_5/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2 
sequential_2/gru_5/zeros/mul/yИ
sequential_2/gru_5/zeros/mulMul)sequential_2/gru_5/strided_slice:output:0'sequential_2/gru_5/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_5/zeros/mul
sequential_2/gru_5/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2!
sequential_2/gru_5/zeros/Less/yГ
sequential_2/gru_5/zeros/LessLess sequential_2/gru_5/zeros/mul:z:0(sequential_2/gru_5/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_5/zeros/Less
!sequential_2/gru_5/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2#
!sequential_2/gru_5/zeros/packed/1Я
sequential_2/gru_5/zeros/packedPack)sequential_2/gru_5/strided_slice:output:0*sequential_2/gru_5/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2!
sequential_2/gru_5/zeros/packed
sequential_2/gru_5/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
sequential_2/gru_5/zeros/ConstС
sequential_2/gru_5/zerosFill(sequential_2/gru_5/zeros/packed:output:0'sequential_2/gru_5/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sequential_2/gru_5/zeros
!sequential_2/gru_5/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2#
!sequential_2/gru_5/transpose/permа
sequential_2/gru_5/transpose	Transpose"sequential_2/gru_4/transpose_1:y:0*sequential_2/gru_5/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
sequential_2/gru_5/transpose
sequential_2/gru_5/Shape_1Shape sequential_2/gru_5/transpose:y:0*
T0*
_output_shapes
:2
sequential_2/gru_5/Shape_1
(sequential_2/gru_5/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(sequential_2/gru_5/strided_slice_1/stackЂ
*sequential_2/gru_5/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_5/strided_slice_1/stack_1Ђ
*sequential_2/gru_5/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_5/strided_slice_1/stack_2р
"sequential_2/gru_5/strided_slice_1StridedSlice#sequential_2/gru_5/Shape_1:output:01sequential_2/gru_5/strided_slice_1/stack:output:03sequential_2/gru_5/strided_slice_1/stack_1:output:03sequential_2/gru_5/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2$
"sequential_2/gru_5/strided_slice_1Ћ
.sequential_2/gru_5/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ20
.sequential_2/gru_5/TensorArrayV2/element_shapeў
 sequential_2/gru_5/TensorArrayV2TensorListReserve7sequential_2/gru_5/TensorArrayV2/element_shape:output:0+sequential_2/gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02"
 sequential_2/gru_5/TensorArrayV2х
Hsequential_2/gru_5/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2J
Hsequential_2/gru_5/TensorArrayUnstack/TensorListFromTensor/element_shapeФ
:sequential_2/gru_5/TensorArrayUnstack/TensorListFromTensorTensorListFromTensor sequential_2/gru_5/transpose:y:0Qsequential_2/gru_5/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02<
:sequential_2/gru_5/TensorArrayUnstack/TensorListFromTensor
(sequential_2/gru_5/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2*
(sequential_2/gru_5/strided_slice_2/stackЂ
*sequential_2/gru_5/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_5/strided_slice_2/stack_1Ђ
*sequential_2/gru_5/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_5/strided_slice_2/stack_2ю
"sequential_2/gru_5/strided_slice_2StridedSlice sequential_2/gru_5/transpose:y:01sequential_2/gru_5/strided_slice_2/stack:output:03sequential_2/gru_5/strided_slice_2/stack_1:output:03sequential_2/gru_5/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2$
"sequential_2/gru_5/strided_slice_2Й
-sequential_2/gru_5/gru_cell_5/ones_like/ShapeShape+sequential_2/gru_5/strided_slice_2:output:0*
T0*
_output_shapes
:2/
-sequential_2/gru_5/gru_cell_5/ones_like/ShapeЃ
-sequential_2/gru_5/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2/
-sequential_2/gru_5/gru_cell_5/ones_like/Constќ
'sequential_2/gru_5/gru_cell_5/ones_likeFill6sequential_2/gru_5/gru_cell_5/ones_like/Shape:output:06sequential_2/gru_5/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/gru_cell_5/ones_likeв
,sequential_2/gru_5/gru_cell_5/ReadVariableOpReadVariableOp5sequential_2_gru_5_gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02.
,sequential_2/gru_5/gru_cell_5/ReadVariableOpФ
%sequential_2/gru_5/gru_cell_5/unstackUnpack4sequential_2/gru_5/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2'
%sequential_2/gru_5/gru_cell_5/unstackо
!sequential_2/gru_5/gru_cell_5/mulMul+sequential_2/gru_5/strided_slice_2:output:00sequential_2/gru_5/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_5/gru_cell_5/mulч
3sequential_2/gru_5/gru_cell_5/MatMul/ReadVariableOpReadVariableOp<sequential_2_gru_5_gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype025
3sequential_2/gru_5/gru_cell_5/MatMul/ReadVariableOpь
$sequential_2/gru_5/gru_cell_5/MatMulMatMul%sequential_2/gru_5/gru_cell_5/mul:z:0;sequential_2/gru_5/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02&
$sequential_2/gru_5/gru_cell_5/MatMulы
%sequential_2/gru_5/gru_cell_5/BiasAddBiasAdd.sequential_2/gru_5/gru_cell_5/MatMul:product:0.sequential_2/gru_5/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02'
%sequential_2/gru_5/gru_cell_5/BiasAdd
#sequential_2/gru_5/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2%
#sequential_2/gru_5/gru_cell_5/ConstЉ
-sequential_2/gru_5/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2/
-sequential_2/gru_5/gru_cell_5/split/split_dimЄ
#sequential_2/gru_5/gru_cell_5/splitSplit6sequential_2/gru_5/gru_cell_5/split/split_dim:output:0.sequential_2/gru_5/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2%
#sequential_2/gru_5/gru_cell_5/splitэ
5sequential_2/gru_5/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp>sequential_2_gru_5_gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype027
5sequential_2/gru_5/gru_cell_5/MatMul_1/ReadVariableOpю
&sequential_2/gru_5/gru_cell_5/MatMul_1MatMul!sequential_2/gru_5/zeros:output:0=sequential_2/gru_5/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02(
&sequential_2/gru_5/gru_cell_5/MatMul_1ё
'sequential_2/gru_5/gru_cell_5/BiasAdd_1BiasAdd0sequential_2/gru_5/gru_cell_5/MatMul_1:product:0.sequential_2/gru_5/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02)
'sequential_2/gru_5/gru_cell_5/BiasAdd_1Ѓ
%sequential_2/gru_5/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2'
%sequential_2/gru_5/gru_cell_5/Const_1­
/sequential_2/gru_5/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ21
/sequential_2/gru_5/gru_cell_5/split_1/split_dimщ
%sequential_2/gru_5/gru_cell_5/split_1SplitV0sequential_2/gru_5/gru_cell_5/BiasAdd_1:output:0.sequential_2/gru_5/gru_cell_5/Const_1:output:08sequential_2/gru_5/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2'
%sequential_2/gru_5/gru_cell_5/split_1п
!sequential_2/gru_5/gru_cell_5/addAddV2,sequential_2/gru_5/gru_cell_5/split:output:0.sequential_2/gru_5/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_5/gru_cell_5/addВ
%sequential_2/gru_5/gru_cell_5/SigmoidSigmoid%sequential_2/gru_5/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%sequential_2/gru_5/gru_cell_5/Sigmoidу
#sequential_2/gru_5/gru_cell_5/add_1AddV2,sequential_2/gru_5/gru_cell_5/split:output:1.sequential_2/gru_5/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/add_1И
'sequential_2/gru_5/gru_cell_5/Sigmoid_1Sigmoid'sequential_2/gru_5/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/gru_cell_5/Sigmoid_1р
#sequential_2/gru_5/gru_cell_5/mul_1Mul+sequential_2/gru_5/gru_cell_5/Sigmoid_1:y:0.sequential_2/gru_5/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/mul_1м
#sequential_2/gru_5/gru_cell_5/add_2AddV2,sequential_2/gru_5/gru_cell_5/split:output:2'sequential_2/gru_5/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/add_2И
'sequential_2/gru_5/gru_cell_5/Sigmoid_2Sigmoid'sequential_2/gru_5/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/gru_cell_5/Sigmoid_2б
#sequential_2/gru_5/gru_cell_5/mul_2Mul)sequential_2/gru_5/gru_cell_5/Sigmoid:y:0!sequential_2/gru_5/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/mul_2
#sequential_2/gru_5/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2%
#sequential_2/gru_5/gru_cell_5/sub/xи
!sequential_2/gru_5/gru_cell_5/subSub,sequential_2/gru_5/gru_cell_5/sub/x:output:0)sequential_2/gru_5/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!sequential_2/gru_5/gru_cell_5/subз
#sequential_2/gru_5/gru_cell_5/mul_3Mul%sequential_2/gru_5/gru_cell_5/sub:z:0+sequential_2/gru_5/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/mul_3з
#sequential_2/gru_5/gru_cell_5/add_3AddV2'sequential_2/gru_5/gru_cell_5/mul_2:z:0'sequential_2/gru_5/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/gru_cell_5/add_3Е
0sequential_2/gru_5/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0sequential_2/gru_5/TensorArrayV2_1/element_shape
"sequential_2/gru_5/TensorArrayV2_1TensorListReserve9sequential_2/gru_5/TensorArrayV2_1/element_shape:output:0+sequential_2/gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02$
"sequential_2/gru_5/TensorArrayV2_1t
sequential_2/gru_5/timeConst*
_output_shapes
: *
dtype0*
value	B : 2
sequential_2/gru_5/timeЅ
+sequential_2/gru_5/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2-
+sequential_2/gru_5/while/maximum_iterations
%sequential_2/gru_5/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2'
%sequential_2/gru_5/while/loop_counterЏ
sequential_2/gru_5/whileWhile.sequential_2/gru_5/while/loop_counter:output:04sequential_2/gru_5/while/maximum_iterations:output:0 sequential_2/gru_5/time:output:0+sequential_2/gru_5/TensorArrayV2_1:handle:0!sequential_2/gru_5/zeros:output:0+sequential_2/gru_5/strided_slice_1:output:0Jsequential_2/gru_5/TensorArrayUnstack/TensorListFromTensor:output_handle:05sequential_2_gru_5_gru_cell_5_readvariableop_resource<sequential_2_gru_5_gru_cell_5_matmul_readvariableop_resource>sequential_2_gru_5_gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*/
body'R%
#sequential_2_gru_5_while_body_67927*/
cond'R%
#sequential_2_gru_5_while_cond_67926*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
sequential_2/gru_5/whileл
Csequential_2/gru_5/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2E
Csequential_2/gru_5/TensorArrayV2Stack/TensorListStack/element_shapeЕ
5sequential_2/gru_5/TensorArrayV2Stack/TensorListStackTensorListStack!sequential_2/gru_5/while:output:3Lsequential_2/gru_5/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype027
5sequential_2/gru_5/TensorArrayV2Stack/TensorListStackЇ
(sequential_2/gru_5/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2*
(sequential_2/gru_5/strided_slice_3/stackЂ
*sequential_2/gru_5/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2,
*sequential_2/gru_5/strided_slice_3/stack_1Ђ
*sequential_2/gru_5/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2,
*sequential_2/gru_5/strided_slice_3/stack_2
"sequential_2/gru_5/strided_slice_3StridedSlice>sequential_2/gru_5/TensorArrayV2Stack/TensorListStack:tensor:01sequential_2/gru_5/strided_slice_3/stack:output:03sequential_2/gru_5/strided_slice_3/stack_1:output:03sequential_2/gru_5/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2$
"sequential_2/gru_5/strided_slice_3
#sequential_2/gru_5/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2%
#sequential_2/gru_5/transpose_1/permђ
sequential_2/gru_5/transpose_1	Transpose>sequential_2/gru_5/TensorArrayV2Stack/TensorListStack:tensor:0,sequential_2/gru_5/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2 
sequential_2/gru_5/transpose_1
sequential_2/gru_5/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
sequential_2/gru_5/runtimeЬ
*sequential_2/dense_2/MatMul/ReadVariableOpReadVariableOp3sequential_2_dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02,
*sequential_2/dense_2/MatMul/ReadVariableOpз
sequential_2/dense_2/MatMulMatMul+sequential_2/gru_5/strided_slice_3:output:02sequential_2/dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sequential_2/dense_2/MatMulЫ
+sequential_2/dense_2/BiasAdd/ReadVariableOpReadVariableOp4sequential_2_dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_2/dense_2/BiasAdd/ReadVariableOpе
sequential_2/dense_2/BiasAddBiasAdd%sequential_2/dense_2/MatMul:product:03sequential_2/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sequential_2/dense_2/BiasAdd 
sequential_2/dense_2/SigmoidSigmoid%sequential_2/dense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sequential_2/dense_2/SigmoidЊ
IdentityIdentity sequential_2/dense_2/Sigmoid:y:0^sequential_2/gru_4/while^sequential_2/gru_5/while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::24
sequential_2/gru_4/whilesequential_2/gru_4/while24
sequential_2/gru_5/whilesequential_2/gru_5/while:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input
Ы
Ѕ
while_cond_69142
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69142___redundant_placeholder03
/while_while_cond_69142___redundant_placeholder13
/while_while_cond_69142___redundant_placeholder23
/while_while_cond_69142___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
њ

gru_4_while_cond_70352(
$gru_4_while_gru_4_while_loop_counter.
*gru_4_while_gru_4_while_maximum_iterations
gru_4_while_placeholder
gru_4_while_placeholder_1
gru_4_while_placeholder_2*
&gru_4_while_less_gru_4_strided_slice_1?
;gru_4_while_gru_4_while_cond_70352___redundant_placeholder0?
;gru_4_while_gru_4_while_cond_70352___redundant_placeholder1?
;gru_4_while_gru_4_while_cond_70352___redundant_placeholder2?
;gru_4_while_gru_4_while_cond_70352___redundant_placeholder3
gru_4_while_identity

gru_4/while/LessLessgru_4_while_placeholder&gru_4_while_less_gru_4_strided_slice_1*
T0*
_output_shapes
: 2
gru_4/while/Lesso
gru_4/while/IdentityIdentitygru_4/while/Less:z:0*
T0
*
_output_shapes
: 2
gru_4/while/Identity"5
gru_4_while_identitygru_4/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
П

G__inference_sequential_2_layer_call_and_return_conditional_losses_70118
embedding_2_input
embedding_2_70095
gru_4_70098
gru_4_70100
gru_4_70102
gru_5_70105
gru_5_70107
gru_5_70109
dense_2_70112
dense_2_70114
identityЂdense_2/StatefulPartitionedCallЂ#embedding_2/StatefulPartitionedCallЂgru_4/StatefulPartitionedCallЂgru_5/StatefulPartitionedCall
#embedding_2/StatefulPartitionedCallStatefulPartitionedCallembedding_2_inputembedding_2_70095*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*#
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *O
fJRH
F__inference_embedding_2_layer_call_and_return_conditional_losses_692302%
#embedding_2/StatefulPartitionedCallМ
gru_4/StatefulPartitionedCallStatefulPartitionedCall,embedding_2/StatefulPartitionedCall:output:0gru_4_70098gru_4_70100gru_4_70102*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_696232
gru_4/StatefulPartitionedCallБ
gru_5/StatefulPartitionedCallStatefulPartitionedCall&gru_4/StatefulPartitionedCall:output:0gru_5_70105gru_5_70107gru_5_70109*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_700342
gru_5/StatefulPartitionedCallЌ
dense_2/StatefulPartitionedCallStatefulPartitionedCall&gru_5/StatefulPartitionedCall:output:0dense_2_70112dense_2_70114*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_700752!
dense_2/StatefulPartitionedCall
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall$^embedding_2/StatefulPartitionedCall^gru_4/StatefulPartitionedCall^gru_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2J
#embedding_2/StatefulPartitionedCall#embedding_2/StatefulPartitionedCall2>
gru_4/StatefulPartitionedCallgru_4/StatefulPartitionedCall2>
gru_5/StatefulPartitionedCallgru_5/StatefulPartitionedCall:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input
Џ

%__inference_gru_4_layer_call_fn_71892
inputs_0
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputs_0unknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_684952
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
И

F__inference_embedding_2_layer_call_and_return_conditional_losses_69230

inputs
embedding_lookup_69224
identity^
CastCastinputs*

DstT0*

SrcT0*(
_output_shapes
:џџџџџџџџџш2
CastЮ
embedding_lookupResourceGatherembedding_lookup_69224Cast:y:0*
Tindices0*)
_class
loc:@embedding_lookup/69224*,
_output_shapes
:џџџџџџџџџш*
dtype02
embedding_lookupП
embedding_lookup/IdentityIdentityembedding_lookup:output:0*
T0*)
_class
loc:@embedding_lookup/69224*,
_output_shapes
:џџџџџџџџџш2
embedding_lookup/IdentityЁ
embedding_lookup/Identity_1Identity"embedding_lookup/Identity:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
embedding_lookup/Identity_1}
IdentityIdentity$embedding_lookup/Identity_1:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*+
_input_shapes
:џџџџџџџџџш::P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ы
Ѕ
while_cond_69939
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69939___redundant_placeholder03
/while_while_cond_69939___redundant_placeholder13
/while_while_cond_69939___redundant_placeholder23
/while_while_cond_69939___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
ы[
я
@__inference_gru_4_layer_call_and_return_conditional_losses_69623

inputs&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_like
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69529*
condR
while_cond_69528*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimep
IdentityIdentitytranspose_1:y:0^while*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Іv
а
gru_5_while_body_70564(
$gru_5_while_gru_5_while_loop_counter.
*gru_5_while_gru_5_while_maximum_iterations
gru_5_while_placeholder
gru_5_while_placeholder_1
gru_5_while_placeholder_2'
#gru_5_while_gru_5_strided_slice_1_0c
_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_04
0gru_5_while_gru_cell_5_readvariableop_resource_0;
7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0=
9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0
gru_5_while_identity
gru_5_while_identity_1
gru_5_while_identity_2
gru_5_while_identity_3
gru_5_while_identity_4%
!gru_5_while_gru_5_strided_slice_1a
]gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor2
.gru_5_while_gru_cell_5_readvariableop_resource9
5gru_5_while_gru_cell_5_matmul_readvariableop_resource;
7gru_5_while_gru_cell_5_matmul_1_readvariableop_resourceЯ
=gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2?
=gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeї
/gru_5/while/TensorArrayV2Read/TensorListGetItemTensorListGetItem_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_0gru_5_while_placeholderFgru_5/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype021
/gru_5/while/TensorArrayV2Read/TensorListGetItemЖ
&gru_5/while/gru_cell_5/ones_like/ShapeShape6gru_5/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2(
&gru_5/while/gru_cell_5/ones_like/Shape
&gru_5/while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2(
&gru_5/while/gru_cell_5/ones_like/Constр
 gru_5/while/gru_cell_5/ones_likeFill/gru_5/while/gru_cell_5/ones_like/Shape:output:0/gru_5/while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/ones_like
$gru_5/while/gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2&
$gru_5/while/gru_cell_5/dropout/Constл
"gru_5/while/gru_cell_5/dropout/MulMul)gru_5/while/gru_cell_5/ones_like:output:0-gru_5/while/gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2$
"gru_5/while/gru_cell_5/dropout/MulЅ
$gru_5/while/gru_cell_5/dropout/ShapeShape)gru_5/while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2&
$gru_5/while/gru_cell_5/dropout/Shape
;gru_5/while/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform-gru_5/while/gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ФЛ2=
;gru_5/while/gru_cell_5/dropout/random_uniform/RandomUniformЃ
-gru_5/while/gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2/
-gru_5/while/gru_cell_5/dropout/GreaterEqual/y
+gru_5/while/gru_cell_5/dropout/GreaterEqualGreaterEqualDgru_5/while/gru_cell_5/dropout/random_uniform/RandomUniform:output:06gru_5/while/gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+gru_5/while/gru_cell_5/dropout/GreaterEqualФ
#gru_5/while/gru_cell_5/dropout/CastCast/gru_5/while/gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2%
#gru_5/while/gru_cell_5/dropout/Castж
$gru_5/while/gru_cell_5/dropout/Mul_1Mul&gru_5/while/gru_cell_5/dropout/Mul:z:0'gru_5/while/gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_5/while/gru_cell_5/dropout/Mul_1
&gru_5/while/gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2(
&gru_5/while/gru_cell_5/dropout_1/Constс
$gru_5/while/gru_cell_5/dropout_1/MulMul)gru_5/while/gru_cell_5/ones_like:output:0/gru_5/while/gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_5/while/gru_cell_5/dropout_1/MulЉ
&gru_5/while/gru_cell_5/dropout_1/ShapeShape)gru_5/while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2(
&gru_5/while/gru_cell_5/dropout_1/Shape
=gru_5/while/gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform/gru_5/while/gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2уЖТ2?
=gru_5/while/gru_cell_5/dropout_1/random_uniform/RandomUniformЇ
/gru_5/while/gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?21
/gru_5/while/gru_cell_5/dropout_1/GreaterEqual/yЂ
-gru_5/while/gru_cell_5/dropout_1/GreaterEqualGreaterEqualFgru_5/while/gru_cell_5/dropout_1/random_uniform/RandomUniform:output:08gru_5/while/gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-gru_5/while/gru_cell_5/dropout_1/GreaterEqualЪ
%gru_5/while/gru_cell_5/dropout_1/CastCast1gru_5/while/gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2'
%gru_5/while/gru_cell_5/dropout_1/Castо
&gru_5/while/gru_cell_5/dropout_1/Mul_1Mul(gru_5/while/gru_cell_5/dropout_1/Mul:z:0)gru_5/while/gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&gru_5/while/gru_cell_5/dropout_1/Mul_1
&gru_5/while/gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2(
&gru_5/while/gru_cell_5/dropout_2/Constс
$gru_5/while/gru_cell_5/dropout_2/MulMul)gru_5/while/gru_cell_5/ones_like:output:0/gru_5/while/gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_5/while/gru_cell_5/dropout_2/MulЉ
&gru_5/while/gru_cell_5/dropout_2/ShapeShape)gru_5/while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2(
&gru_5/while/gru_cell_5/dropout_2/Shape
=gru_5/while/gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform/gru_5/while/gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Ђѕ2?
=gru_5/while/gru_cell_5/dropout_2/random_uniform/RandomUniformЇ
/gru_5/while/gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?21
/gru_5/while/gru_cell_5/dropout_2/GreaterEqual/yЂ
-gru_5/while/gru_cell_5/dropout_2/GreaterEqualGreaterEqualFgru_5/while/gru_cell_5/dropout_2/random_uniform/RandomUniform:output:08gru_5/while/gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-gru_5/while/gru_cell_5/dropout_2/GreaterEqualЪ
%gru_5/while/gru_cell_5/dropout_2/CastCast1gru_5/while/gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2'
%gru_5/while/gru_cell_5/dropout_2/Castо
&gru_5/while/gru_cell_5/dropout_2/Mul_1Mul(gru_5/while/gru_cell_5/dropout_2/Mul:z:0)gru_5/while/gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&gru_5/while/gru_cell_5/dropout_2/Mul_1П
%gru_5/while/gru_cell_5/ReadVariableOpReadVariableOp0gru_5_while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02'
%gru_5/while/gru_cell_5/ReadVariableOpЏ
gru_5/while/gru_cell_5/unstackUnpack-gru_5/while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2 
gru_5/while/gru_cell_5/unstackг
gru_5/while/gru_cell_5/mulMul6gru_5/while/TensorArrayV2Read/TensorListGetItem:item:0(gru_5/while/gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mulд
,gru_5/while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02.
,gru_5/while/gru_cell_5/MatMul/ReadVariableOpа
gru_5/while/gru_cell_5/MatMulMatMulgru_5/while/gru_cell_5/mul:z:04gru_5/while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/while/gru_cell_5/MatMulЯ
gru_5/while/gru_cell_5/BiasAddBiasAdd'gru_5/while/gru_cell_5/MatMul:product:0'gru_5/while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02 
gru_5/while/gru_cell_5/BiasAdd~
gru_5/while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/gru_cell_5/Const
&gru_5/while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2(
&gru_5/while/gru_cell_5/split/split_dim
gru_5/while/gru_cell_5/splitSplit/gru_5/while/gru_cell_5/split/split_dim:output:0'gru_5/while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/while/gru_cell_5/splitк
.gru_5/while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype020
.gru_5/while/gru_cell_5/MatMul_1/ReadVariableOpб
gru_5/while/gru_cell_5/MatMul_1MatMulgru_5_while_placeholder_26gru_5/while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02!
gru_5/while/gru_cell_5/MatMul_1е
 gru_5/while/gru_cell_5/BiasAdd_1BiasAdd)gru_5/while/gru_cell_5/MatMul_1:product:0'gru_5/while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02"
 gru_5/while/gru_cell_5/BiasAdd_1
gru_5/while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2 
gru_5/while/gru_cell_5/Const_1
(gru_5/while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2*
(gru_5/while/gru_cell_5/split_1/split_dimЦ
gru_5/while/gru_cell_5/split_1SplitV)gru_5/while/gru_cell_5/BiasAdd_1:output:0'gru_5/while/gru_cell_5/Const_1:output:01gru_5/while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2 
gru_5/while/gru_cell_5/split_1У
gru_5/while/gru_cell_5/addAddV2%gru_5/while/gru_cell_5/split:output:0'gru_5/while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add
gru_5/while/gru_cell_5/SigmoidSigmoidgru_5/while/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_5/while/gru_cell_5/SigmoidЧ
gru_5/while/gru_cell_5/add_1AddV2%gru_5/while/gru_cell_5/split:output:1'gru_5/while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_1Ѓ
 gru_5/while/gru_cell_5/Sigmoid_1Sigmoid gru_5/while/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/Sigmoid_1Ф
gru_5/while/gru_cell_5/mul_1Mul$gru_5/while/gru_cell_5/Sigmoid_1:y:0'gru_5/while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_1Р
gru_5/while/gru_cell_5/add_2AddV2%gru_5/while/gru_cell_5/split:output:2 gru_5/while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_2Ѓ
 gru_5/while/gru_cell_5/Sigmoid_2Sigmoid gru_5/while/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/Sigmoid_2Д
gru_5/while/gru_cell_5/mul_2Mul"gru_5/while/gru_cell_5/Sigmoid:y:0gru_5_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_2
gru_5/while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_5/while/gru_cell_5/sub/xМ
gru_5/while/gru_cell_5/subSub%gru_5/while/gru_cell_5/sub/x:output:0"gru_5/while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/subЛ
gru_5/while/gru_cell_5/mul_3Mulgru_5/while/gru_cell_5/sub:z:0$gru_5/while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_3Л
gru_5/while/gru_cell_5/add_3AddV2 gru_5/while/gru_cell_5/mul_2:z:0 gru_5/while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_3ќ
0gru_5/while/TensorArrayV2Write/TensorListSetItemTensorListSetItemgru_5_while_placeholder_1gru_5_while_placeholder gru_5/while/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype022
0gru_5/while/TensorArrayV2Write/TensorListSetItemh
gru_5/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/add/y
gru_5/while/addAddV2gru_5_while_placeholdergru_5/while/add/y:output:0*
T0*
_output_shapes
: 2
gru_5/while/addl
gru_5/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/add_1/y
gru_5/while/add_1AddV2$gru_5_while_gru_5_while_loop_countergru_5/while/add_1/y:output:0*
T0*
_output_shapes
: 2
gru_5/while/add_1p
gru_5/while/IdentityIdentitygru_5/while/add_1:z:0*
T0*
_output_shapes
: 2
gru_5/while/Identity
gru_5/while/Identity_1Identity*gru_5_while_gru_5_while_maximum_iterations*
T0*
_output_shapes
: 2
gru_5/while/Identity_1r
gru_5/while/Identity_2Identitygru_5/while/add:z:0*
T0*
_output_shapes
: 2
gru_5/while/Identity_2
gru_5/while/Identity_3Identity@gru_5/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
gru_5/while/Identity_3
gru_5/while/Identity_4Identity gru_5/while/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/Identity_4"H
!gru_5_while_gru_5_strided_slice_1#gru_5_while_gru_5_strided_slice_1_0"t
7gru_5_while_gru_cell_5_matmul_1_readvariableop_resource9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0"p
5gru_5_while_gru_cell_5_matmul_readvariableop_resource7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0"b
.gru_5_while_gru_cell_5_readvariableop_resource0gru_5_while_gru_cell_5_readvariableop_resource_0"5
gru_5_while_identitygru_5/while/Identity:output:0"9
gru_5_while_identity_1gru_5/while/Identity_1:output:0"9
gru_5_while_identity_2gru_5/while/Identity_2:output:0"9
gru_5_while_identity_3gru_5/while/Identity_3:output:0"9
gru_5_while_identity_4gru_5/while/Identity_4:output:0"Р
]gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 


G__inference_sequential_2_layer_call_and_return_conditional_losses_70196

inputs
embedding_2_70173
gru_4_70176
gru_4_70178
gru_4_70180
gru_5_70183
gru_5_70185
gru_5_70187
dense_2_70190
dense_2_70192
identityЂdense_2/StatefulPartitionedCallЂ#embedding_2/StatefulPartitionedCallЂgru_4/StatefulPartitionedCallЂgru_5/StatefulPartitionedCall
#embedding_2/StatefulPartitionedCallStatefulPartitionedCallinputsembedding_2_70173*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*#
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *O
fJRH
F__inference_embedding_2_layer_call_and_return_conditional_losses_692302%
#embedding_2/StatefulPartitionedCallМ
gru_4/StatefulPartitionedCallStatefulPartitionedCall,embedding_2/StatefulPartitionedCall:output:0gru_4_70176gru_4_70178gru_4_70180*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_696232
gru_4/StatefulPartitionedCallБ
gru_5/StatefulPartitionedCallStatefulPartitionedCall&gru_4/StatefulPartitionedCall:output:0gru_5_70183gru_5_70185gru_5_70187*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_700342
gru_5/StatefulPartitionedCallЌ
dense_2/StatefulPartitionedCallStatefulPartitionedCall&gru_5/StatefulPartitionedCall:output:0dense_2_70190dense_2_70192*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_700752!
dense_2/StatefulPartitionedCall
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall$^embedding_2/StatefulPartitionedCall^gru_4/StatefulPartitionedCall^gru_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2J
#embedding_2/StatefulPartitionedCall#embedding_2/StatefulPartitionedCall2>
gru_4/StatefulPartitionedCallgru_4/StatefulPartitionedCall2>
gru_5/StatefulPartitionedCallgru_5/StatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
И

F__inference_embedding_2_layer_call_and_return_conditional_losses_71088

inputs
embedding_lookup_71082
identity^
CastCastinputs*

DstT0*

SrcT0*(
_output_shapes
:џџџџџџџџџш2
CastЮ
embedding_lookupResourceGatherembedding_lookup_71082Cast:y:0*
Tindices0*)
_class
loc:@embedding_lookup/71082*,
_output_shapes
:џџџџџџџџџш*
dtype02
embedding_lookupП
embedding_lookup/IdentityIdentityembedding_lookup:output:0*
T0*)
_class
loc:@embedding_lookup/71082*,
_output_shapes
:џџџџџџџџџш2
embedding_lookup/IdentityЁ
embedding_lookup/Identity_1Identity"embedding_lookup/Identity:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
embedding_lookup/Identity_1}
IdentityIdentity$embedding_lookup/Identity_1:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*+
_input_shapes
:џџџџџџџџџш::P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ь<
Ю
@__inference_gru_4_layer_call_and_return_conditional_losses_68613

inputs
gru_cell_4_68537
gru_cell_4_68539
gru_cell_4_68541
identityЂ"gru_cell_4/StatefulPartitionedCallЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputstranspose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2ц
"gru_cell_4/StatefulPartitionedCallStatefulPartitionedCallstrided_slice_2:output:0zeros:output:0gru_cell_4_68537gru_cell_4_68539gru_cell_4_68541*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681722$
"gru_cell_4/StatefulPartitionedCall
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterп
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0gru_cell_4_68537gru_cell_4_68539gru_cell_4_68541*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_68549*
condR
while_cond_68548*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtime
IdentityIdentitytranspose_1:y:0#^gru_cell_4/StatefulPartitionedCall^while*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2H
"gru_cell_4/StatefulPartitionedCall"gru_cell_4/StatefulPartitionedCall2
whilewhile:\ X
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
вi
Ў
while_body_71596
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like
while/gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_4/dropout/ConstУ
while/gru_cell_4/dropout/MulMul#while/gru_cell_4/ones_like:output:0'while/gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/Mul
while/gru_cell_4/dropout/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_4/dropout/Shape
5while/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Умi27
5while/gru_cell_4/dropout/random_uniform/RandomUniform
'while/gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_4/dropout/GreaterEqual/y
%while/gru_cell_4/dropout/GreaterEqualGreaterEqual>while/gru_cell_4/dropout/random_uniform/RandomUniform:output:00while/gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_4/dropout/GreaterEqualВ
while/gru_cell_4/dropout/CastCast)while/gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/CastО
while/gru_cell_4/dropout/Mul_1Mul while/gru_cell_4/dropout/Mul:z:0!while/gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout/Mul_1
 while/gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_1/ConstЩ
while/gru_cell_4/dropout_1/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_1/Mul
 while/gru_cell_4/dropout_1/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_1/Shape
7while/gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ц§ъ29
7while/gru_cell_4/dropout_1/random_uniform/RandomUniform
)while/gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_1/GreaterEqual/y
'while/gru_cell_4/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_1/GreaterEqualИ
while/gru_cell_4/dropout_1/CastCast+while/gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_1/CastЦ
 while/gru_cell_4/dropout_1/Mul_1Mul"while/gru_cell_4/dropout_1/Mul:z:0#while/gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_1/Mul_1
 while/gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_2/ConstЩ
while/gru_cell_4/dropout_2/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_2/Mul
 while/gru_cell_4/dropout_2/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_2/Shape
7while/gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2­Ю29
7while/gru_cell_4/dropout_2/random_uniform/RandomUniform
)while/gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_2/GreaterEqual/y
'while/gru_cell_4/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_2/GreaterEqualИ
while/gru_cell_4/dropout_2/CastCast+while/gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_2/CastЦ
 while/gru_cell_4/dropout_2/Mul_1Mul"while/gru_cell_4/dropout_2/Mul:z:0#while/gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_2/Mul_1­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackЛ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
вi
Ў
while_body_69749
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like
while/gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_5/dropout/ConstУ
while/gru_cell_5/dropout/MulMul#while/gru_cell_5/ones_like:output:0'while/gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/Mul
while/gru_cell_5/dropout/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_5/dropout/Shape
5while/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2С27
5while/gru_cell_5/dropout/random_uniform/RandomUniform
'while/gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_5/dropout/GreaterEqual/y
%while/gru_cell_5/dropout/GreaterEqualGreaterEqual>while/gru_cell_5/dropout/random_uniform/RandomUniform:output:00while/gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_5/dropout/GreaterEqualВ
while/gru_cell_5/dropout/CastCast)while/gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/CastО
while/gru_cell_5/dropout/Mul_1Mul while/gru_cell_5/dropout/Mul:z:0!while/gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout/Mul_1
 while/gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_1/ConstЩ
while/gru_cell_5/dropout_1/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_1/Mul
 while/gru_cell_5/dropout_1/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_1/Shape
7while/gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Ј29
7while/gru_cell_5/dropout_1/random_uniform/RandomUniform
)while/gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_1/GreaterEqual/y
'while/gru_cell_5/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_1/GreaterEqualИ
while/gru_cell_5/dropout_1/CastCast+while/gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_1/CastЦ
 while/gru_cell_5/dropout_1/Mul_1Mul"while/gru_cell_5/dropout_1/Mul:z:0#while/gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_1/Mul_1
 while/gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_2/ConstЩ
while/gru_cell_5/dropout_2/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_2/Mul
 while/gru_cell_5/dropout_2/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_2/Shape
7while/gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Є229
7while/gru_cell_5/dropout_2/random_uniform/RandomUniform
)while/gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_2/GreaterEqual/y
'while/gru_cell_5/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_2/GreaterEqualИ
while/gru_cell_5/dropout_2/CastCast+while/gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_2/CastЦ
 while/gru_cell_5/dropout_2/Mul_1Mul"while/gru_cell_5/dropout_2/Mul:z:0#while/gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_2/Mul_1­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackЛ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Љ
Њ
B__inference_dense_2_layer_call_and_return_conditional_losses_70075

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoid_
IdentityIdentitySigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ:::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
к	
Ќ
*__inference_gru_cell_4_layer_call_fn_72871

inputs
states_0
unknown
	unknown_0
	unknown_1
identity

identity_1ЂStatefulPartitionedCallЁ
StatefulPartitionedCallStatefulPartitionedCallinputsstates_0unknown	unknown_0	unknown_1*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681722
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
Ы
Ѕ
while_cond_71191
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_71191___redundant_placeholder03
/while_while_cond_71191___redundant_placeholder13
/while_while_cond_71191___redundant_placeholder23
/while_while_cond_71191___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
њ

gru_4_while_cond_70767(
$gru_4_while_gru_4_while_loop_counter.
*gru_4_while_gru_4_while_maximum_iterations
gru_4_while_placeholder
gru_4_while_placeholder_1
gru_4_while_placeholder_2*
&gru_4_while_less_gru_4_strided_slice_1?
;gru_4_while_gru_4_while_cond_70767___redundant_placeholder0?
;gru_4_while_gru_4_while_cond_70767___redundant_placeholder1?
;gru_4_while_gru_4_while_cond_70767___redundant_placeholder2?
;gru_4_while_gru_4_while_cond_70767___redundant_placeholder3
gru_4_while_identity

gru_4/while/LessLessgru_4_while_placeholder&gru_4_while_less_gru_4_strided_slice_1*
T0*
_output_shapes
: 2
gru_4/while/Lesso
gru_4/while/IdentityIdentitygru_4/while/Less:z:0*
T0
*
_output_shapes
: 2
gru_4/while/Identity"5
gru_4_while_identitygru_4/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
ц
ѕ
,__inference_sequential_2_layer_call_fn_70217
embedding_2_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
identityЂStatefulPartitionedCallн
StatefulPartitionedCallStatefulPartitionedCallembedding_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_701962
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input
Ы
Ѕ
while_cond_71382
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_71382___redundant_placeholder03
/while_while_cond_71382___redundant_placeholder13
/while_while_cond_71382___redundant_placeholder23
/while_while_cond_71382___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Ш<
Ю
@__inference_gru_5_layer_call_and_return_conditional_losses_69207

inputs
gru_cell_5_69131
gru_cell_5_69133
gru_cell_5_69135
identityЂ"gru_cell_5/StatefulPartitionedCallЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputstranspose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2ц
"gru_cell_5/StatefulPartitionedCallStatefulPartitionedCallstrided_slice_2:output:0zeros:output:0gru_cell_5_69131gru_cell_5_69133gru_cell_5_69135*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687662$
"gru_cell_5/StatefulPartitionedCall
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterп
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0gru_cell_5_69131gru_cell_5_69133gru_cell_5_69135*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69143*
condR
while_cond_69142*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtime
IdentityIdentitystrided_slice_3:output:0#^gru_cell_5/StatefulPartitionedCall^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2H
"gru_cell_5/StatefulPartitionedCall"gru_cell_5/StatefulPartitionedCall2
whilewhile:\ X
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
}
я
@__inference_gru_4_layer_call_and_return_conditional_losses_69456

inputs&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_likey
gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout/ConstЋ
gru_cell_4/dropout/MulMulgru_cell_4/ones_like:output:0!gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul
gru_cell_4/dropout/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout/Shapeє
/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2квс21
/gru_cell_4/dropout/random_uniform/RandomUniform
!gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_4/dropout/GreaterEqual/yъ
gru_cell_4/dropout/GreaterEqualGreaterEqual8gru_cell_4/dropout/random_uniform/RandomUniform:output:0*gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_4/dropout/GreaterEqual 
gru_cell_4/dropout/CastCast#gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/CastІ
gru_cell_4/dropout/Mul_1Mulgru_cell_4/dropout/Mul:z:0gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul_1}
gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_1/ConstБ
gru_cell_4/dropout_1/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul
gru_cell_4/dropout_1/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_1/Shapeњ
1gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2цЌю23
1gru_cell_4/dropout_1/random_uniform/RandomUniform
#gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_1/GreaterEqual/yђ
!gru_cell_4/dropout_1/GreaterEqualGreaterEqual:gru_cell_4/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_1/GreaterEqualІ
gru_cell_4/dropout_1/CastCast%gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/CastЎ
gru_cell_4/dropout_1/Mul_1Mulgru_cell_4/dropout_1/Mul:z:0gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul_1}
gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_2/ConstБ
gru_cell_4/dropout_2/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul
gru_cell_4/dropout_2/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_2/Shapeњ
1gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2­ п23
1gru_cell_4/dropout_2/random_uniform/RandomUniform
#gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_2/GreaterEqual/yђ
!gru_cell_4/dropout_2/GreaterEqualGreaterEqual:gru_cell_4/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_2/GreaterEqualІ
gru_cell_4/dropout_2/CastCast%gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/CastЎ
gru_cell_4/dropout_2/Mul_1Mulgru_cell_4/dropout_2/Mul:z:0gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul_1
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69338*
condR
while_cond_69337*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimep
IdentityIdentitytranspose_1:y:0^while*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ы
Ѕ
while_cond_72594
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_72594___redundant_placeholder03
/while_while_cond_72594___redundant_placeholder13
/while_while_cond_72594___redundant_placeholder23
/while_while_cond_72594___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
 \
ё
@__inference_gru_5_layer_call_and_return_conditional_losses_72285
inputs_0&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileF
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputs_0transpose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_like
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_72191*
condR
while_cond_72190*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2
whilewhile:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
џ

%__inference_gru_5_layer_call_fn_72711

inputs
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall§
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_700342
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
к	
Ќ
*__inference_gru_cell_4_layer_call_fn_72857

inputs
states_0
unknown
	unknown_0
	unknown_1
identity

identity_1ЂStatefulPartitionedCallЁ
StatefulPartitionedCallStatefulPartitionedCallinputsstates_0unknown	unknown_0	unknown_1*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681282
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
ы[
я
@__inference_gru_4_layer_call_and_return_conditional_losses_71477

inputs&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_like
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_71383*
condR
while_cond_71382*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimep
IdentityIdentitytranspose_1:y:0^while*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ы
Ѕ
while_cond_69337
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69337___redundant_placeholder03
/while_while_cond_69337___redundant_placeholder13
/while_while_cond_69337___redundant_placeholder23
/while_while_cond_69337___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
є

#sequential_2_gru_4_while_cond_67763B
>sequential_2_gru_4_while_sequential_2_gru_4_while_loop_counterH
Dsequential_2_gru_4_while_sequential_2_gru_4_while_maximum_iterations(
$sequential_2_gru_4_while_placeholder*
&sequential_2_gru_4_while_placeholder_1*
&sequential_2_gru_4_while_placeholder_2D
@sequential_2_gru_4_while_less_sequential_2_gru_4_strided_slice_1Y
Usequential_2_gru_4_while_sequential_2_gru_4_while_cond_67763___redundant_placeholder0Y
Usequential_2_gru_4_while_sequential_2_gru_4_while_cond_67763___redundant_placeholder1Y
Usequential_2_gru_4_while_sequential_2_gru_4_while_cond_67763___redundant_placeholder2Y
Usequential_2_gru_4_while_sequential_2_gru_4_while_cond_67763___redundant_placeholder3%
!sequential_2_gru_4_while_identity
Я
sequential_2/gru_4/while/LessLess$sequential_2_gru_4_while_placeholder@sequential_2_gru_4_while_less_sequential_2_gru_4_strided_slice_1*
T0*
_output_shapes
: 2
sequential_2/gru_4/while/Less
!sequential_2/gru_4/while/IdentityIdentity!sequential_2/gru_4/while/Less:z:0*
T0
*
_output_shapes
: 2#
!sequential_2/gru_4/while/Identity"O
!sequential_2_gru_4_while_identity*sequential_2/gru_4/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Ы
Ѕ
while_cond_71595
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_71595___redundant_placeholder03
/while_while_cond_71595___redundant_placeholder13
/while_while_cond_71595___redundant_placeholder23
/while_while_cond_71595___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:


G__inference_sequential_2_layer_call_and_return_conditional_losses_70147

inputs
embedding_2_70124
gru_4_70127
gru_4_70129
gru_4_70131
gru_5_70134
gru_5_70136
gru_5_70138
dense_2_70141
dense_2_70143
identityЂdense_2/StatefulPartitionedCallЂ#embedding_2/StatefulPartitionedCallЂgru_4/StatefulPartitionedCallЂgru_5/StatefulPartitionedCall
#embedding_2/StatefulPartitionedCallStatefulPartitionedCallinputsembedding_2_70124*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*#
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *O
fJRH
F__inference_embedding_2_layer_call_and_return_conditional_losses_692302%
#embedding_2/StatefulPartitionedCallМ
gru_4/StatefulPartitionedCallStatefulPartitionedCall,embedding_2/StatefulPartitionedCall:output:0gru_4_70127gru_4_70129gru_4_70131*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_694562
gru_4/StatefulPartitionedCallБ
gru_5/StatefulPartitionedCallStatefulPartitionedCall&gru_4/StatefulPartitionedCall:output:0gru_5_70134gru_5_70136gru_5_70138*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_698672
gru_5/StatefulPartitionedCallЌ
dense_2/StatefulPartitionedCallStatefulPartitionedCall&gru_5/StatefulPartitionedCall:output:0dense_2_70141dense_2_70143*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_700752!
dense_2/StatefulPartitionedCall
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall$^embedding_2/StatefulPartitionedCall^gru_4/StatefulPartitionedCall^gru_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2J
#embedding_2/StatefulPartitionedCall#embedding_2/StatefulPartitionedCall2>
gru_4/StatefulPartitionedCallgru_4/StatefulPartitionedCall2>
gru_5/StatefulPartitionedCallgru_5/StatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
ђD
Ў
while_body_69940
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackМ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
ђD
Ў
while_body_71383
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackМ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
џ

%__inference_gru_5_layer_call_fn_72700

inputs
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall§
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_698672
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
њ

gru_5_while_cond_70930(
$gru_5_while_gru_5_while_loop_counter.
*gru_5_while_gru_5_while_maximum_iterations
gru_5_while_placeholder
gru_5_while_placeholder_1
gru_5_while_placeholder_2*
&gru_5_while_less_gru_5_strided_slice_1?
;gru_5_while_gru_5_while_cond_70930___redundant_placeholder0?
;gru_5_while_gru_5_while_cond_70930___redundant_placeholder1?
;gru_5_while_gru_5_while_cond_70930___redundant_placeholder2?
;gru_5_while_gru_5_while_cond_70930___redundant_placeholder3
gru_5_while_identity

gru_5/while/LessLessgru_5_while_placeholder&gru_5_while_less_gru_5_strided_slice_1*
T0*
_output_shapes
: 2
gru_5/while/Lesso
gru_5/while/IdentityIdentitygru_5/while/Less:z:0*
T0
*
_output_shapes
: 2
gru_5/while/Identity"5
gru_5_while_identitygru_5/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Џ

%__inference_gru_4_layer_call_fn_71903
inputs_0
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputs_0unknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_686132
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
Ы
Ѕ
while_cond_68548
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_68548___redundant_placeholder03
/while_while_cond_68548___redundant_placeholder13
/while_while_cond_68548___redundant_placeholder23
/while_while_cond_68548___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Є\
ё
@__inference_gru_4_layer_call_and_return_conditional_losses_71881
inputs_0&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileF
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputs_0transpose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_like
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_71787*
condR
while_cond_71786*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimex
IdentityIdentitytranspose_1:y:0^while*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2
whilewhile:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
гi
Ў
while_body_72000
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like
while/gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_5/dropout/ConstУ
while/gru_cell_5/dropout/MulMul#while/gru_cell_5/ones_like:output:0'while/gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/Mul
while/gru_cell_5/dropout/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_5/dropout/Shape
5while/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2С27
5while/gru_cell_5/dropout/random_uniform/RandomUniform
'while/gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_5/dropout/GreaterEqual/y
%while/gru_cell_5/dropout/GreaterEqualGreaterEqual>while/gru_cell_5/dropout/random_uniform/RandomUniform:output:00while/gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_5/dropout/GreaterEqualВ
while/gru_cell_5/dropout/CastCast)while/gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/dropout/CastО
while/gru_cell_5/dropout/Mul_1Mul while/gru_cell_5/dropout/Mul:z:0!while/gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout/Mul_1
 while/gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_1/ConstЩ
while/gru_cell_5/dropout_1/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_1/Mul
 while/gru_cell_5/dropout_1/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_1/Shape
7while/gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЈГЈ29
7while/gru_cell_5/dropout_1/random_uniform/RandomUniform
)while/gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_1/GreaterEqual/y
'while/gru_cell_5/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_1/GreaterEqualИ
while/gru_cell_5/dropout_1/CastCast+while/gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_1/CastЦ
 while/gru_cell_5/dropout_1/Mul_1Mul"while/gru_cell_5/dropout_1/Mul:z:0#while/gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_1/Mul_1
 while/gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_5/dropout_2/ConstЩ
while/gru_cell_5/dropout_2/MulMul#while/gru_cell_5/ones_like:output:0)while/gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_5/dropout_2/Mul
 while/gru_cell_5/dropout_2/ShapeShape#while/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/dropout_2/Shape
7while/gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2њ 29
7while/gru_cell_5/dropout_2/random_uniform/RandomUniform
)while/gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_5/dropout_2/GreaterEqual/y
'while/gru_cell_5/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_5/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_5/dropout_2/GreaterEqualИ
while/gru_cell_5/dropout_2/CastCast+while/gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_5/dropout_2/CastЦ
 while/gru_cell_5/dropout_2/Mul_1Mul"while/gru_cell_5/dropout_2/Mul:z:0#while/gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_5/dropout_2/Mul_1­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackЛ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
и
|
'__inference_dense_2_layer_call_fn_72731

inputs
unknown
	unknown_0
identityЂStatefulPartitionedCallђ
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_700752
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs


!__inference__traced_restore_73260
file_prefix+
'assignvariableop_embedding_2_embeddings%
!assignvariableop_1_dense_2_kernel#
assignvariableop_2_dense_2_bias 
assignvariableop_3_adam_iter"
assignvariableop_4_adam_beta_1"
assignvariableop_5_adam_beta_2!
assignvariableop_6_adam_decay)
%assignvariableop_7_adam_learning_rate.
*assignvariableop_8_gru_4_gru_cell_4_kernel8
4assignvariableop_9_gru_4_gru_cell_4_recurrent_kernel-
)assignvariableop_10_gru_4_gru_cell_4_bias/
+assignvariableop_11_gru_5_gru_cell_5_kernel9
5assignvariableop_12_gru_5_gru_cell_5_recurrent_kernel-
)assignvariableop_13_gru_5_gru_cell_5_bias
assignvariableop_14_total
assignvariableop_15_count
assignvariableop_16_total_1
assignvariableop_17_count_15
1assignvariableop_18_adam_embedding_2_embeddings_m-
)assignvariableop_19_adam_dense_2_kernel_m+
'assignvariableop_20_adam_dense_2_bias_m6
2assignvariableop_21_adam_gru_4_gru_cell_4_kernel_m@
<assignvariableop_22_adam_gru_4_gru_cell_4_recurrent_kernel_m4
0assignvariableop_23_adam_gru_4_gru_cell_4_bias_m6
2assignvariableop_24_adam_gru_5_gru_cell_5_kernel_m@
<assignvariableop_25_adam_gru_5_gru_cell_5_recurrent_kernel_m4
0assignvariableop_26_adam_gru_5_gru_cell_5_bias_m5
1assignvariableop_27_adam_embedding_2_embeddings_v-
)assignvariableop_28_adam_dense_2_kernel_v+
'assignvariableop_29_adam_dense_2_bias_v6
2assignvariableop_30_adam_gru_4_gru_cell_4_kernel_v@
<assignvariableop_31_adam_gru_4_gru_cell_4_recurrent_kernel_v4
0assignvariableop_32_adam_gru_4_gru_cell_4_bias_v6
2assignvariableop_33_adam_gru_5_gru_cell_5_kernel_v@
<assignvariableop_34_adam_gru_5_gru_cell_5_recurrent_kernel_v4
0assignvariableop_35_adam_gru_5_gru_cell_5_bias_v
identity_37ЂAssignVariableOpЂAssignVariableOp_1ЂAssignVariableOp_10ЂAssignVariableOp_11ЂAssignVariableOp_12ЂAssignVariableOp_13ЂAssignVariableOp_14ЂAssignVariableOp_15ЂAssignVariableOp_16ЂAssignVariableOp_17ЂAssignVariableOp_18ЂAssignVariableOp_19ЂAssignVariableOp_2ЂAssignVariableOp_20ЂAssignVariableOp_21ЂAssignVariableOp_22ЂAssignVariableOp_23ЂAssignVariableOp_24ЂAssignVariableOp_25ЂAssignVariableOp_26ЂAssignVariableOp_27ЂAssignVariableOp_28ЂAssignVariableOp_29ЂAssignVariableOp_3ЂAssignVariableOp_30ЂAssignVariableOp_31ЂAssignVariableOp_32ЂAssignVariableOp_33ЂAssignVariableOp_34ЂAssignVariableOp_35ЂAssignVariableOp_4ЂAssignVariableOp_5ЂAssignVariableOp_6ЂAssignVariableOp_7ЂAssignVariableOp_8ЂAssignVariableOp_9о
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:%*
dtype0*ъ
valueрBн%B:layer_with_weights-0/embeddings/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_namesи
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:%*
dtype0*]
valueTBR%B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slicesч
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*Њ
_output_shapes
:::::::::::::::::::::::::::::::::::::*3
dtypes)
'2%	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

IdentityІ
AssignVariableOpAssignVariableOp'assignvariableop_embedding_2_embeddingsIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1І
AssignVariableOp_1AssignVariableOp!assignvariableop_1_dense_2_kernelIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2Є
AssignVariableOp_2AssignVariableOpassignvariableop_2_dense_2_biasIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0	*
_output_shapes
:2

Identity_3Ё
AssignVariableOp_3AssignVariableOpassignvariableop_3_adam_iterIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype0	2
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4Ѓ
AssignVariableOp_4AssignVariableOpassignvariableop_4_adam_beta_1Identity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5Ѓ
AssignVariableOp_5AssignVariableOpassignvariableop_5_adam_beta_2Identity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6Ђ
AssignVariableOp_6AssignVariableOpassignvariableop_6_adam_decayIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7Њ
AssignVariableOp_7AssignVariableOp%assignvariableop_7_adam_learning_rateIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8Џ
AssignVariableOp_8AssignVariableOp*assignvariableop_8_gru_4_gru_cell_4_kernelIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9Й
AssignVariableOp_9AssignVariableOp4assignvariableop_9_gru_4_gru_cell_4_recurrent_kernelIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10Б
AssignVariableOp_10AssignVariableOp)assignvariableop_10_gru_4_gru_cell_4_biasIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11Г
AssignVariableOp_11AssignVariableOp+assignvariableop_11_gru_5_gru_cell_5_kernelIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12Н
AssignVariableOp_12AssignVariableOp5assignvariableop_12_gru_5_gru_cell_5_recurrent_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13Б
AssignVariableOp_13AssignVariableOp)assignvariableop_13_gru_5_gru_cell_5_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14Ё
AssignVariableOp_14AssignVariableOpassignvariableop_14_totalIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15Ё
AssignVariableOp_15AssignVariableOpassignvariableop_15_countIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16Ѓ
AssignVariableOp_16AssignVariableOpassignvariableop_16_total_1Identity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17Ѓ
AssignVariableOp_17AssignVariableOpassignvariableop_17_count_1Identity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18Й
AssignVariableOp_18AssignVariableOp1assignvariableop_18_adam_embedding_2_embeddings_mIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19Б
AssignVariableOp_19AssignVariableOp)assignvariableop_19_adam_dense_2_kernel_mIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_19n
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:2
Identity_20Џ
AssignVariableOp_20AssignVariableOp'assignvariableop_20_adam_dense_2_bias_mIdentity_20:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_20n
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:2
Identity_21К
AssignVariableOp_21AssignVariableOp2assignvariableop_21_adam_gru_4_gru_cell_4_kernel_mIdentity_21:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_21n
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:2
Identity_22Ф
AssignVariableOp_22AssignVariableOp<assignvariableop_22_adam_gru_4_gru_cell_4_recurrent_kernel_mIdentity_22:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_22n
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:2
Identity_23И
AssignVariableOp_23AssignVariableOp0assignvariableop_23_adam_gru_4_gru_cell_4_bias_mIdentity_23:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_23n
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:2
Identity_24К
AssignVariableOp_24AssignVariableOp2assignvariableop_24_adam_gru_5_gru_cell_5_kernel_mIdentity_24:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_24n
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:2
Identity_25Ф
AssignVariableOp_25AssignVariableOp<assignvariableop_25_adam_gru_5_gru_cell_5_recurrent_kernel_mIdentity_25:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_25n
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:2
Identity_26И
AssignVariableOp_26AssignVariableOp0assignvariableop_26_adam_gru_5_gru_cell_5_bias_mIdentity_26:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_26n
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:2
Identity_27Й
AssignVariableOp_27AssignVariableOp1assignvariableop_27_adam_embedding_2_embeddings_vIdentity_27:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_27n
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:2
Identity_28Б
AssignVariableOp_28AssignVariableOp)assignvariableop_28_adam_dense_2_kernel_vIdentity_28:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_28n
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:2
Identity_29Џ
AssignVariableOp_29AssignVariableOp'assignvariableop_29_adam_dense_2_bias_vIdentity_29:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_29n
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:2
Identity_30К
AssignVariableOp_30AssignVariableOp2assignvariableop_30_adam_gru_4_gru_cell_4_kernel_vIdentity_30:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_30n
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:2
Identity_31Ф
AssignVariableOp_31AssignVariableOp<assignvariableop_31_adam_gru_4_gru_cell_4_recurrent_kernel_vIdentity_31:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_31n
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:2
Identity_32И
AssignVariableOp_32AssignVariableOp0assignvariableop_32_adam_gru_4_gru_cell_4_bias_vIdentity_32:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_32n
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:2
Identity_33К
AssignVariableOp_33AssignVariableOp2assignvariableop_33_adam_gru_5_gru_cell_5_kernel_vIdentity_33:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_33n
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:2
Identity_34Ф
AssignVariableOp_34AssignVariableOp<assignvariableop_34_adam_gru_5_gru_cell_5_recurrent_kernel_vIdentity_34:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_34n
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:2
Identity_35И
AssignVariableOp_35AssignVariableOp0assignvariableop_35_adam_gru_5_gru_cell_5_bias_vIdentity_35:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_359
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOpі
Identity_36Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_36щ
Identity_37IdentityIdentity_36:output:0^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*
T0*
_output_shapes
: 2
Identity_37"#
identity_37Identity_37:output:0*Ї
_input_shapes
: ::::::::::::::::::::::::::::::::::::2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
Ы
Ѕ
while_cond_69748
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69748___redundant_placeholder03
/while_while_cond_69748___redundant_placeholder13
/while_while_cond_69748___redundant_placeholder23
/while_while_cond_69748___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
;
ч
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_68722

inputs

states
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likec
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Const
dropout/MulMulones_like:output:0dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul`
dropout/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout/Shapeг
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЁЎЬ2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/yО
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul_1g
dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_1/Const
dropout_1/MulMulones_like:output:0dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Muld
dropout_1/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_1/Shapeй
&dropout_1/random_uniform/RandomUniformRandomUniformdropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2дсЬ2(
&dropout_1/random_uniform/RandomUniformy
dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_1/GreaterEqual/yЦ
dropout_1/GreaterEqualGreaterEqual/dropout_1/random_uniform/RandomUniform:output:0!dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/GreaterEqual
dropout_1/CastCastdropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Cast
dropout_1/Mul_1Muldropout_1/Mul:z:0dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Mul_1g
dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_2/Const
dropout_2/MulMulones_like:output:0dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Muld
dropout_2/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_2/Shapeй
&dropout_2/random_uniform/RandomUniformRandomUniformdropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЗнА2(
&dropout_2/random_uniform/RandomUniformy
dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_2/GreaterEqual/yЦ
dropout_2/GreaterEqualGreaterEqual/dropout_2/random_uniform/RandomUniform:output:0!dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/GreaterEqual
dropout_2/CastCastdropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Cast
dropout_2/Mul_1Muldropout_2/Mul:z:0dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Mul_1x
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack^
mulMulinputsdropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOpy
MatMul_1MatMulstatesMatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2\
mul_2MulSigmoid:y:0states*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:OK
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_namestates
Ы
Ѕ
while_cond_69024
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69024___redundant_placeholder03
/while_while_cond_69024___redundant_placeholder13
/while_while_cond_69024___redundant_placeholder23
/while_while_cond_69024___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Ь
q
+__inference_embedding_2_layer_call_fn_71095

inputs
unknown
identityЂStatefulPartitionedCallю
StatefulPartitionedCallStatefulPartitionedCallinputsunknown*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*#
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *O
fJRH
F__inference_embedding_2_layer_call_and_return_conditional_losses_692302
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*+
_input_shapes
:џџџџџџџџџш:22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ѕv
а
gru_4_while_body_70353(
$gru_4_while_gru_4_while_loop_counter.
*gru_4_while_gru_4_while_maximum_iterations
gru_4_while_placeholder
gru_4_while_placeholder_1
gru_4_while_placeholder_2'
#gru_4_while_gru_4_strided_slice_1_0c
_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_04
0gru_4_while_gru_cell_4_readvariableop_resource_0;
7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0=
9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0
gru_4_while_identity
gru_4_while_identity_1
gru_4_while_identity_2
gru_4_while_identity_3
gru_4_while_identity_4%
!gru_4_while_gru_4_strided_slice_1a
]gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor2
.gru_4_while_gru_cell_4_readvariableop_resource9
5gru_4_while_gru_cell_4_matmul_readvariableop_resource;
7gru_4_while_gru_cell_4_matmul_1_readvariableop_resourceЯ
=gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2?
=gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeї
/gru_4/while/TensorArrayV2Read/TensorListGetItemTensorListGetItem_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_0gru_4_while_placeholderFgru_4/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype021
/gru_4/while/TensorArrayV2Read/TensorListGetItemЖ
&gru_4/while/gru_cell_4/ones_like/ShapeShape6gru_4/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2(
&gru_4/while/gru_cell_4/ones_like/Shape
&gru_4/while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2(
&gru_4/while/gru_cell_4/ones_like/Constр
 gru_4/while/gru_cell_4/ones_likeFill/gru_4/while/gru_cell_4/ones_like/Shape:output:0/gru_4/while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/ones_like
$gru_4/while/gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2&
$gru_4/while/gru_cell_4/dropout/Constл
"gru_4/while/gru_cell_4/dropout/MulMul)gru_4/while/gru_cell_4/ones_like:output:0-gru_4/while/gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2$
"gru_4/while/gru_cell_4/dropout/MulЅ
$gru_4/while/gru_cell_4/dropout/ShapeShape)gru_4/while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2&
$gru_4/while/gru_cell_4/dropout/Shape
;gru_4/while/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform-gru_4/while/gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2"2=
;gru_4/while/gru_cell_4/dropout/random_uniform/RandomUniformЃ
-gru_4/while/gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2/
-gru_4/while/gru_cell_4/dropout/GreaterEqual/y
+gru_4/while/gru_cell_4/dropout/GreaterEqualGreaterEqualDgru_4/while/gru_cell_4/dropout/random_uniform/RandomUniform:output:06gru_4/while/gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+gru_4/while/gru_cell_4/dropout/GreaterEqualФ
#gru_4/while/gru_cell_4/dropout/CastCast/gru_4/while/gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2%
#gru_4/while/gru_cell_4/dropout/Castж
$gru_4/while/gru_cell_4/dropout/Mul_1Mul&gru_4/while/gru_cell_4/dropout/Mul:z:0'gru_4/while/gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_4/while/gru_cell_4/dropout/Mul_1
&gru_4/while/gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2(
&gru_4/while/gru_cell_4/dropout_1/Constс
$gru_4/while/gru_cell_4/dropout_1/MulMul)gru_4/while/gru_cell_4/ones_like:output:0/gru_4/while/gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_4/while/gru_cell_4/dropout_1/MulЉ
&gru_4/while/gru_cell_4/dropout_1/ShapeShape)gru_4/while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2(
&gru_4/while/gru_cell_4/dropout_1/Shape
=gru_4/while/gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform/gru_4/while/gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЈЎп2?
=gru_4/while/gru_cell_4/dropout_1/random_uniform/RandomUniformЇ
/gru_4/while/gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?21
/gru_4/while/gru_cell_4/dropout_1/GreaterEqual/yЂ
-gru_4/while/gru_cell_4/dropout_1/GreaterEqualGreaterEqualFgru_4/while/gru_cell_4/dropout_1/random_uniform/RandomUniform:output:08gru_4/while/gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-gru_4/while/gru_cell_4/dropout_1/GreaterEqualЪ
%gru_4/while/gru_cell_4/dropout_1/CastCast1gru_4/while/gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2'
%gru_4/while/gru_cell_4/dropout_1/Castо
&gru_4/while/gru_cell_4/dropout_1/Mul_1Mul(gru_4/while/gru_cell_4/dropout_1/Mul:z:0)gru_4/while/gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&gru_4/while/gru_cell_4/dropout_1/Mul_1
&gru_4/while/gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2(
&gru_4/while/gru_cell_4/dropout_2/Constс
$gru_4/while/gru_cell_4/dropout_2/MulMul)gru_4/while/gru_cell_4/ones_like:output:0/gru_4/while/gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2&
$gru_4/while/gru_cell_4/dropout_2/MulЉ
&gru_4/while/gru_cell_4/dropout_2/ShapeShape)gru_4/while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2(
&gru_4/while/gru_cell_4/dropout_2/Shape
=gru_4/while/gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform/gru_4/while/gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2тўЬ2?
=gru_4/while/gru_cell_4/dropout_2/random_uniform/RandomUniformЇ
/gru_4/while/gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?21
/gru_4/while/gru_cell_4/dropout_2/GreaterEqual/yЂ
-gru_4/while/gru_cell_4/dropout_2/GreaterEqualGreaterEqualFgru_4/while/gru_cell_4/dropout_2/random_uniform/RandomUniform:output:08gru_4/while/gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-gru_4/while/gru_cell_4/dropout_2/GreaterEqualЪ
%gru_4/while/gru_cell_4/dropout_2/CastCast1gru_4/while/gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2'
%gru_4/while/gru_cell_4/dropout_2/Castо
&gru_4/while/gru_cell_4/dropout_2/Mul_1Mul(gru_4/while/gru_cell_4/dropout_2/Mul:z:0)gru_4/while/gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2(
&gru_4/while/gru_cell_4/dropout_2/Mul_1П
%gru_4/while/gru_cell_4/ReadVariableOpReadVariableOp0gru_4_while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02'
%gru_4/while/gru_cell_4/ReadVariableOpЏ
gru_4/while/gru_cell_4/unstackUnpack-gru_4/while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2 
gru_4/while/gru_cell_4/unstackг
gru_4/while/gru_cell_4/mulMul6gru_4/while/TensorArrayV2Read/TensorListGetItem:item:0(gru_4/while/gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mulд
,gru_4/while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02.
,gru_4/while/gru_cell_4/MatMul/ReadVariableOpа
gru_4/while/gru_cell_4/MatMulMatMulgru_4/while/gru_cell_4/mul:z:04gru_4/while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/while/gru_cell_4/MatMulЯ
gru_4/while/gru_cell_4/BiasAddBiasAdd'gru_4/while/gru_cell_4/MatMul:product:0'gru_4/while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02 
gru_4/while/gru_cell_4/BiasAdd~
gru_4/while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/gru_cell_4/Const
&gru_4/while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2(
&gru_4/while/gru_cell_4/split/split_dim
gru_4/while/gru_cell_4/splitSplit/gru_4/while/gru_cell_4/split/split_dim:output:0'gru_4/while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/while/gru_cell_4/splitк
.gru_4/while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype020
.gru_4/while/gru_cell_4/MatMul_1/ReadVariableOpб
gru_4/while/gru_cell_4/MatMul_1MatMulgru_4_while_placeholder_26gru_4/while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02!
gru_4/while/gru_cell_4/MatMul_1е
 gru_4/while/gru_cell_4/BiasAdd_1BiasAdd)gru_4/while/gru_cell_4/MatMul_1:product:0'gru_4/while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02"
 gru_4/while/gru_cell_4/BiasAdd_1
gru_4/while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2 
gru_4/while/gru_cell_4/Const_1
(gru_4/while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2*
(gru_4/while/gru_cell_4/split_1/split_dimЦ
gru_4/while/gru_cell_4/split_1SplitV)gru_4/while/gru_cell_4/BiasAdd_1:output:0'gru_4/while/gru_cell_4/Const_1:output:01gru_4/while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2 
gru_4/while/gru_cell_4/split_1У
gru_4/while/gru_cell_4/addAddV2%gru_4/while/gru_cell_4/split:output:0'gru_4/while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add
gru_4/while/gru_cell_4/SigmoidSigmoidgru_4/while/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_4/while/gru_cell_4/SigmoidЧ
gru_4/while/gru_cell_4/add_1AddV2%gru_4/while/gru_cell_4/split:output:1'gru_4/while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_1Ѓ
 gru_4/while/gru_cell_4/Sigmoid_1Sigmoid gru_4/while/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/Sigmoid_1Ф
gru_4/while/gru_cell_4/mul_1Mul$gru_4/while/gru_cell_4/Sigmoid_1:y:0'gru_4/while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_1Р
gru_4/while/gru_cell_4/add_2AddV2%gru_4/while/gru_cell_4/split:output:2 gru_4/while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_2Ѓ
 gru_4/while/gru_cell_4/Sigmoid_2Sigmoid gru_4/while/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/Sigmoid_2Д
gru_4/while/gru_cell_4/mul_2Mul"gru_4/while/gru_cell_4/Sigmoid:y:0gru_4_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_2
gru_4/while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_4/while/gru_cell_4/sub/xМ
gru_4/while/gru_cell_4/subSub%gru_4/while/gru_cell_4/sub/x:output:0"gru_4/while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/subЛ
gru_4/while/gru_cell_4/mul_3Mulgru_4/while/gru_cell_4/sub:z:0$gru_4/while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_3Л
gru_4/while/gru_cell_4/add_3AddV2 gru_4/while/gru_cell_4/mul_2:z:0 gru_4/while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_3ќ
0gru_4/while/TensorArrayV2Write/TensorListSetItemTensorListSetItemgru_4_while_placeholder_1gru_4_while_placeholder gru_4/while/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype022
0gru_4/while/TensorArrayV2Write/TensorListSetItemh
gru_4/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/add/y
gru_4/while/addAddV2gru_4_while_placeholdergru_4/while/add/y:output:0*
T0*
_output_shapes
: 2
gru_4/while/addl
gru_4/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/add_1/y
gru_4/while/add_1AddV2$gru_4_while_gru_4_while_loop_countergru_4/while/add_1/y:output:0*
T0*
_output_shapes
: 2
gru_4/while/add_1p
gru_4/while/IdentityIdentitygru_4/while/add_1:z:0*
T0*
_output_shapes
: 2
gru_4/while/Identity
gru_4/while/Identity_1Identity*gru_4_while_gru_4_while_maximum_iterations*
T0*
_output_shapes
: 2
gru_4/while/Identity_1r
gru_4/while/Identity_2Identitygru_4/while/add:z:0*
T0*
_output_shapes
: 2
gru_4/while/Identity_2
gru_4/while/Identity_3Identity@gru_4/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
gru_4/while/Identity_3
gru_4/while/Identity_4Identity gru_4/while/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/Identity_4"H
!gru_4_while_gru_4_strided_slice_1#gru_4_while_gru_4_strided_slice_1_0"t
7gru_4_while_gru_cell_4_matmul_1_readvariableop_resource9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0"p
5gru_4_while_gru_cell_4_matmul_readvariableop_resource7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0"b
.gru_4_while_gru_cell_4_readvariableop_resource0gru_4_while_gru_cell_4_readvariableop_resource_0"5
gru_4_while_identitygru_4/while/Identity:output:0"9
gru_4_while_identity_1gru_4/while/Identity_1:output:0"9
gru_4_while_identity_2gru_4/while/Identity_2:output:0"9
gru_4_while_identity_3gru_4/while/Identity_3:output:0"9
gru_4_while_identity_4gru_4/while/Identity_4:output:0"Р
]gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
гi
Ў
while_body_69338
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like
while/gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_4/dropout/ConstУ
while/gru_cell_4/dropout/MulMul#while/gru_cell_4/ones_like:output:0'while/gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/Mul
while/gru_cell_4/dropout/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_4/dropout/Shape
5while/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2юОЩ27
5while/gru_cell_4/dropout/random_uniform/RandomUniform
'while/gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_4/dropout/GreaterEqual/y
%while/gru_cell_4/dropout/GreaterEqualGreaterEqual>while/gru_cell_4/dropout/random_uniform/RandomUniform:output:00while/gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_4/dropout/GreaterEqualВ
while/gru_cell_4/dropout/CastCast)while/gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/CastО
while/gru_cell_4/dropout/Mul_1Mul while/gru_cell_4/dropout/Mul:z:0!while/gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout/Mul_1
 while/gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_1/ConstЩ
while/gru_cell_4/dropout_1/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_1/Mul
 while/gru_cell_4/dropout_1/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_1/Shape
7while/gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2МЕЅ29
7while/gru_cell_4/dropout_1/random_uniform/RandomUniform
)while/gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_1/GreaterEqual/y
'while/gru_cell_4/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_1/GreaterEqualИ
while/gru_cell_4/dropout_1/CastCast+while/gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_1/CastЦ
 while/gru_cell_4/dropout_1/Mul_1Mul"while/gru_cell_4/dropout_1/Mul:z:0#while/gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_1/Mul_1
 while/gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_2/ConstЩ
while/gru_cell_4/dropout_2/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_2/Mul
 while/gru_cell_4/dropout_2/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_2/Shape
7while/gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2мо29
7while/gru_cell_4/dropout_2/random_uniform/RandomUniform
)while/gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_2/GreaterEqual/y
'while/gru_cell_4/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_2/GreaterEqualИ
while/gru_cell_4/dropout_2/CastCast+while/gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_2/CastЦ
 while/gru_cell_4/dropout_2/Mul_1Mul"while/gru_cell_4/dropout_2/Mul:z:0#while/gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_2/Mul_1­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackЛ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Ы
Ѕ
while_cond_68430
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_68430___redundant_placeholder03
/while_while_cond_68430___redundant_placeholder13
/while_while_cond_68430___redundant_placeholder23
/while_while_cond_68430___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Л}
ё
@__inference_gru_5_layer_call_and_return_conditional_losses_72118
inputs_0&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileF
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputs_0transpose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_likey
gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout/ConstЋ
gru_cell_5/dropout/MulMulgru_cell_5/ones_like:output:0!gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul
gru_cell_5/dropout/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout/Shapeє
/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2йз21
/gru_cell_5/dropout/random_uniform/RandomUniform
!gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_5/dropout/GreaterEqual/yъ
gru_cell_5/dropout/GreaterEqualGreaterEqual8gru_cell_5/dropout/random_uniform/RandomUniform:output:0*gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_5/dropout/GreaterEqual 
gru_cell_5/dropout/CastCast#gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/CastІ
gru_cell_5/dropout/Mul_1Mulgru_cell_5/dropout/Mul:z:0gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul_1}
gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_1/ConstБ
gru_cell_5/dropout_1/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul
gru_cell_5/dropout_1/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_1/Shapeњ
1gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ѓї23
1gru_cell_5/dropout_1/random_uniform/RandomUniform
#gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_1/GreaterEqual/yђ
!gru_cell_5/dropout_1/GreaterEqualGreaterEqual:gru_cell_5/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_1/GreaterEqualІ
gru_cell_5/dropout_1/CastCast%gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/CastЎ
gru_cell_5/dropout_1/Mul_1Mulgru_cell_5/dropout_1/Mul:z:0gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul_1}
gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_2/ConstБ
gru_cell_5/dropout_2/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul
gru_cell_5/dropout_2/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_2/Shapeљ
1gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2.23
1gru_cell_5/dropout_2/random_uniform/RandomUniform
#gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_2/GreaterEqual/yђ
!gru_cell_5/dropout_2/GreaterEqualGreaterEqual:gru_cell_5/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_2/GreaterEqualІ
gru_cell_5/dropout_2/CastCast%gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/CastЎ
gru_cell_5/dropout_2/Mul_1Mulgru_cell_5/dropout_2/Mul:z:0gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul_1
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_72000*
condR
while_cond_71999*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2
whilewhile:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
}
я
@__inference_gru_5_layer_call_and_return_conditional_losses_72522

inputs&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_likey
gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout/ConstЋ
gru_cell_5/dropout/MulMulgru_cell_5/ones_like:output:0!gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul
gru_cell_5/dropout/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout/Shapeѓ
/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Ќг921
/gru_cell_5/dropout/random_uniform/RandomUniform
!gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_5/dropout/GreaterEqual/yъ
gru_cell_5/dropout/GreaterEqualGreaterEqual8gru_cell_5/dropout/random_uniform/RandomUniform:output:0*gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_5/dropout/GreaterEqual 
gru_cell_5/dropout/CastCast#gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/CastІ
gru_cell_5/dropout/Mul_1Mulgru_cell_5/dropout/Mul:z:0gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul_1}
gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_1/ConstБ
gru_cell_5/dropout_1/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul
gru_cell_5/dropout_1/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_1/Shapeњ
1gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЏПє23
1gru_cell_5/dropout_1/random_uniform/RandomUniform
#gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_1/GreaterEqual/yђ
!gru_cell_5/dropout_1/GreaterEqualGreaterEqual:gru_cell_5/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_1/GreaterEqualІ
gru_cell_5/dropout_1/CastCast%gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/CastЎ
gru_cell_5/dropout_1/Mul_1Mulgru_cell_5/dropout_1/Mul:z:0gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul_1}
gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_2/ConstБ
gru_cell_5/dropout_2/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul
gru_cell_5/dropout_2/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_2/Shapeњ
1gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЪЬў23
1gru_cell_5/dropout_2/random_uniform/RandomUniform
#gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_2/GreaterEqual/yђ
!gru_cell_5/dropout_2/GreaterEqualGreaterEqual:gru_cell_5/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_2/GreaterEqualІ
gru_cell_5/dropout_2/CastCast%gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/CastЎ
gru_cell_5/dropout_2/Mul_1Mulgru_cell_5/dropout_2/Mul:z:0gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul_1
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_72404*
condR
while_cond_72403*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
А 
щ
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72843

inputs
states_0
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likex
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack_
mulMulinputsones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOp{
MatMul_1MatMulstates_0MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2^
mul_2MulSigmoid:y:0states_0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
ўЃ
З
G__inference_sequential_2_layer_call_and_return_conditional_losses_70689

inputs&
"embedding_2_embedding_lookup_70254,
(gru_4_gru_cell_4_readvariableop_resource3
/gru_4_gru_cell_4_matmul_readvariableop_resource5
1gru_4_gru_cell_4_matmul_1_readvariableop_resource,
(gru_5_gru_cell_5_readvariableop_resource3
/gru_5_gru_cell_5_matmul_readvariableop_resource5
1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
&dense_2_matmul_readvariableop_resource+
'dense_2_biasadd_readvariableop_resource
identityЂgru_4/whileЂgru_5/whilev
embedding_2/CastCastinputs*

DstT0*

SrcT0*(
_output_shapes
:џџџџџџџџџш2
embedding_2/Cast
embedding_2/embedding_lookupResourceGather"embedding_2_embedding_lookup_70254embedding_2/Cast:y:0*
Tindices0*5
_class+
)'loc:@embedding_2/embedding_lookup/70254*,
_output_shapes
:џџџџџџџџџш*
dtype02
embedding_2/embedding_lookupя
%embedding_2/embedding_lookup/IdentityIdentity%embedding_2/embedding_lookup:output:0*
T0*5
_class+
)'loc:@embedding_2/embedding_lookup/70254*,
_output_shapes
:џџџџџџџџџш2'
%embedding_2/embedding_lookup/IdentityХ
'embedding_2/embedding_lookup/Identity_1Identity.embedding_2/embedding_lookup/Identity:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2)
'embedding_2/embedding_lookup/Identity_1z
gru_4/ShapeShape0embedding_2/embedding_lookup/Identity_1:output:0*
T0*
_output_shapes
:2
gru_4/Shape
gru_4/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice/stack
gru_4/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice/stack_1
gru_4/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice/stack_2
gru_4/strided_sliceStridedSlicegru_4/Shape:output:0"gru_4/strided_slice/stack:output:0$gru_4/strided_slice/stack_1:output:0$gru_4/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_4/strided_sliceh
gru_4/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/zeros/mul/y
gru_4/zeros/mulMulgru_4/strided_slice:output:0gru_4/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
gru_4/zeros/mulk
gru_4/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
gru_4/zeros/Less/y
gru_4/zeros/LessLessgru_4/zeros/mul:z:0gru_4/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
gru_4/zeros/Lessn
gru_4/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
gru_4/zeros/packed/1
gru_4/zeros/packedPackgru_4/strided_slice:output:0gru_4/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
gru_4/zeros/packedk
gru_4/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_4/zeros/Const
gru_4/zerosFillgru_4/zeros/packed:output:0gru_4/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/zeros
gru_4/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_4/transpose/permЗ
gru_4/transpose	Transpose0embedding_2/embedding_lookup/Identity_1:output:0gru_4/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
gru_4/transposea
gru_4/Shape_1Shapegru_4/transpose:y:0*
T0*
_output_shapes
:2
gru_4/Shape_1
gru_4/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_1/stack
gru_4/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_1/stack_1
gru_4/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_1/stack_2
gru_4/strided_slice_1StridedSlicegru_4/Shape_1:output:0$gru_4/strided_slice_1/stack:output:0&gru_4/strided_slice_1/stack_1:output:0&gru_4/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_4/strided_slice_1
!gru_4/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2#
!gru_4/TensorArrayV2/element_shapeЪ
gru_4/TensorArrayV2TensorListReserve*gru_4/TensorArrayV2/element_shape:output:0gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_4/TensorArrayV2Ы
;gru_4/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2=
;gru_4/TensorArrayUnstack/TensorListFromTensor/element_shape
-gru_4/TensorArrayUnstack/TensorListFromTensorTensorListFromTensorgru_4/transpose:y:0Dgru_4/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02/
-gru_4/TensorArrayUnstack/TensorListFromTensor
gru_4/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_2/stack
gru_4/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_2/stack_1
gru_4/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_2/stack_2 
gru_4/strided_slice_2StridedSlicegru_4/transpose:y:0$gru_4/strided_slice_2/stack:output:0&gru_4/strided_slice_2/stack_1:output:0&gru_4/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_4/strided_slice_2
 gru_4/gru_cell_4/ones_like/ShapeShapegru_4/strided_slice_2:output:0*
T0*
_output_shapes
:2"
 gru_4/gru_cell_4/ones_like/Shape
 gru_4/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 gru_4/gru_cell_4/ones_like/ConstШ
gru_4/gru_cell_4/ones_likeFill)gru_4/gru_cell_4/ones_like/Shape:output:0)gru_4/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/ones_like
gru_4/gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
gru_4/gru_cell_4/dropout/ConstУ
gru_4/gru_cell_4/dropout/MulMul#gru_4/gru_cell_4/ones_like:output:0'gru_4/gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/dropout/Mul
gru_4/gru_cell_4/dropout/ShapeShape#gru_4/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2 
gru_4/gru_cell_4/dropout/Shape
5gru_4/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform'gru_4/gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2уЋ27
5gru_4/gru_cell_4/dropout/random_uniform/RandomUniform
'gru_4/gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'gru_4/gru_cell_4/dropout/GreaterEqual/y
%gru_4/gru_cell_4/dropout/GreaterEqualGreaterEqual>gru_4/gru_cell_4/dropout/random_uniform/RandomUniform:output:00gru_4/gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%gru_4/gru_cell_4/dropout/GreaterEqualВ
gru_4/gru_cell_4/dropout/CastCast)gru_4/gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/dropout/CastО
gru_4/gru_cell_4/dropout/Mul_1Mul gru_4/gru_cell_4/dropout/Mul:z:0!gru_4/gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_4/gru_cell_4/dropout/Mul_1
 gru_4/gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 gru_4/gru_cell_4/dropout_1/ConstЩ
gru_4/gru_cell_4/dropout_1/MulMul#gru_4/gru_cell_4/ones_like:output:0)gru_4/gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_4/gru_cell_4/dropout_1/Mul
 gru_4/gru_cell_4/dropout_1/ShapeShape#gru_4/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 gru_4/gru_cell_4/dropout_1/Shape
7gru_4/gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform)gru_4/gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2џЎж29
7gru_4/gru_cell_4/dropout_1/random_uniform/RandomUniform
)gru_4/gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)gru_4/gru_cell_4/dropout_1/GreaterEqual/y
'gru_4/gru_cell_4/dropout_1/GreaterEqualGreaterEqual@gru_4/gru_cell_4/dropout_1/random_uniform/RandomUniform:output:02gru_4/gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'gru_4/gru_cell_4/dropout_1/GreaterEqualИ
gru_4/gru_cell_4/dropout_1/CastCast+gru_4/gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
gru_4/gru_cell_4/dropout_1/CastЦ
 gru_4/gru_cell_4/dropout_1/Mul_1Mul"gru_4/gru_cell_4/dropout_1/Mul:z:0#gru_4/gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/gru_cell_4/dropout_1/Mul_1
 gru_4/gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 gru_4/gru_cell_4/dropout_2/ConstЩ
gru_4/gru_cell_4/dropout_2/MulMul#gru_4/gru_cell_4/ones_like:output:0)gru_4/gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_4/gru_cell_4/dropout_2/Mul
 gru_4/gru_cell_4/dropout_2/ShapeShape#gru_4/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 gru_4/gru_cell_4/dropout_2/Shape
7gru_4/gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform)gru_4/gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ѕп29
7gru_4/gru_cell_4/dropout_2/random_uniform/RandomUniform
)gru_4/gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)gru_4/gru_cell_4/dropout_2/GreaterEqual/y
'gru_4/gru_cell_4/dropout_2/GreaterEqualGreaterEqual@gru_4/gru_cell_4/dropout_2/random_uniform/RandomUniform:output:02gru_4/gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'gru_4/gru_cell_4/dropout_2/GreaterEqualИ
gru_4/gru_cell_4/dropout_2/CastCast+gru_4/gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
gru_4/gru_cell_4/dropout_2/CastЦ
 gru_4/gru_cell_4/dropout_2/Mul_1Mul"gru_4/gru_cell_4/dropout_2/Mul:z:0#gru_4/gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/gru_cell_4/dropout_2/Mul_1Ћ
gru_4/gru_cell_4/ReadVariableOpReadVariableOp(gru_4_gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02!
gru_4/gru_cell_4/ReadVariableOp
gru_4/gru_cell_4/unstackUnpack'gru_4/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_4/gru_cell_4/unstackЉ
gru_4/gru_cell_4/mulMulgru_4/strided_slice_2:output:0"gru_4/gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mulР
&gru_4/gru_cell_4/MatMul/ReadVariableOpReadVariableOp/gru_4_gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02(
&gru_4/gru_cell_4/MatMul/ReadVariableOpИ
gru_4/gru_cell_4/MatMulMatMulgru_4/gru_cell_4/mul:z:0.gru_4/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/MatMulЗ
gru_4/gru_cell_4/BiasAddBiasAdd!gru_4/gru_cell_4/MatMul:product:0!gru_4/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/BiasAddr
gru_4/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/gru_cell_4/Const
 gru_4/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 gru_4/gru_cell_4/split/split_dim№
gru_4/gru_cell_4/splitSplit)gru_4/gru_cell_4/split/split_dim:output:0!gru_4/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/gru_cell_4/splitЦ
(gru_4/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp1gru_4_gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02*
(gru_4/gru_cell_4/MatMul_1/ReadVariableOpК
gru_4/gru_cell_4/MatMul_1MatMulgru_4/zeros:output:00gru_4/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/MatMul_1Н
gru_4/gru_cell_4/BiasAdd_1BiasAdd#gru_4/gru_cell_4/MatMul_1:product:0!gru_4/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/gru_cell_4/BiasAdd_1
gru_4/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_4/gru_cell_4/Const_1
"gru_4/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"gru_4/gru_cell_4/split_1/split_dimЈ
gru_4/gru_cell_4/split_1SplitV#gru_4/gru_cell_4/BiasAdd_1:output:0!gru_4/gru_cell_4/Const_1:output:0+gru_4/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/gru_cell_4/split_1Ћ
gru_4/gru_cell_4/addAddV2gru_4/gru_cell_4/split:output:0!gru_4/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add
gru_4/gru_cell_4/SigmoidSigmoidgru_4/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/SigmoidЏ
gru_4/gru_cell_4/add_1AddV2gru_4/gru_cell_4/split:output:1!gru_4/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_1
gru_4/gru_cell_4/Sigmoid_1Sigmoidgru_4/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/Sigmoid_1Ќ
gru_4/gru_cell_4/mul_1Mulgru_4/gru_cell_4/Sigmoid_1:y:0!gru_4/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_1Ј
gru_4/gru_cell_4/add_2AddV2gru_4/gru_cell_4/split:output:2gru_4/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_2
gru_4/gru_cell_4/Sigmoid_2Sigmoidgru_4/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/Sigmoid_2
gru_4/gru_cell_4/mul_2Mulgru_4/gru_cell_4/Sigmoid:y:0gru_4/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_2u
gru_4/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_4/gru_cell_4/sub/xЄ
gru_4/gru_cell_4/subSubgru_4/gru_cell_4/sub/x:output:0gru_4/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/subЃ
gru_4/gru_cell_4/mul_3Mulgru_4/gru_cell_4/sub:z:0gru_4/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/mul_3Ѓ
gru_4/gru_cell_4/add_3AddV2gru_4/gru_cell_4/mul_2:z:0gru_4/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/gru_cell_4/add_3
#gru_4/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2%
#gru_4/TensorArrayV2_1/element_shapeа
gru_4/TensorArrayV2_1TensorListReserve,gru_4/TensorArrayV2_1/element_shape:output:0gru_4/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_4/TensorArrayV2_1Z

gru_4/timeConst*
_output_shapes
: *
dtype0*
value	B : 2

gru_4/time
gru_4/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2 
gru_4/while/maximum_iterationsv
gru_4/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
gru_4/while/loop_counterљ
gru_4/whileWhile!gru_4/while/loop_counter:output:0'gru_4/while/maximum_iterations:output:0gru_4/time:output:0gru_4/TensorArrayV2_1:handle:0gru_4/zeros:output:0gru_4/strided_slice_1:output:0=gru_4/TensorArrayUnstack/TensorListFromTensor:output_handle:0(gru_4_gru_cell_4_readvariableop_resource/gru_4_gru_cell_4_matmul_readvariableop_resource1gru_4_gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*"
bodyR
gru_4_while_body_70353*"
condR
gru_4_while_cond_70352*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
gru_4/whileС
6gru_4/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   28
6gru_4/TensorArrayV2Stack/TensorListStack/element_shape
(gru_4/TensorArrayV2Stack/TensorListStackTensorListStackgru_4/while:output:3?gru_4/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02*
(gru_4/TensorArrayV2Stack/TensorListStack
gru_4/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
gru_4/strided_slice_3/stack
gru_4/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
gru_4/strided_slice_3/stack_1
gru_4/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_4/strided_slice_3/stack_2О
gru_4/strided_slice_3StridedSlice1gru_4/TensorArrayV2Stack/TensorListStack:tensor:0$gru_4/strided_slice_3/stack:output:0&gru_4/strided_slice_3/stack_1:output:0&gru_4/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_4/strided_slice_3
gru_4/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_4/transpose_1/permО
gru_4/transpose_1	Transpose1gru_4/TensorArrayV2Stack/TensorListStack:tensor:0gru_4/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
gru_4/transpose_1r
gru_4/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_4/runtime_
gru_5/ShapeShapegru_4/transpose_1:y:0*
T0*
_output_shapes
:2
gru_5/Shape
gru_5/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice/stack
gru_5/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice/stack_1
gru_5/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice/stack_2
gru_5/strided_sliceStridedSlicegru_5/Shape:output:0"gru_5/strided_slice/stack:output:0$gru_5/strided_slice/stack_1:output:0$gru_5/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_5/strided_sliceh
gru_5/zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/zeros/mul/y
gru_5/zeros/mulMulgru_5/strided_slice:output:0gru_5/zeros/mul/y:output:0*
T0*
_output_shapes
: 2
gru_5/zeros/mulk
gru_5/zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
gru_5/zeros/Less/y
gru_5/zeros/LessLessgru_5/zeros/mul:z:0gru_5/zeros/Less/y:output:0*
T0*
_output_shapes
: 2
gru_5/zeros/Lessn
gru_5/zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
gru_5/zeros/packed/1
gru_5/zeros/packedPackgru_5/strided_slice:output:0gru_5/zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
gru_5/zeros/packedk
gru_5/zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_5/zeros/Const
gru_5/zerosFillgru_5/zeros/packed:output:0gru_5/zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/zeros
gru_5/transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_5/transpose/perm
gru_5/transpose	Transposegru_4/transpose_1:y:0gru_5/transpose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
gru_5/transposea
gru_5/Shape_1Shapegru_5/transpose:y:0*
T0*
_output_shapes
:2
gru_5/Shape_1
gru_5/strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_1/stack
gru_5/strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_1/stack_1
gru_5/strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_1/stack_2
gru_5/strided_slice_1StridedSlicegru_5/Shape_1:output:0$gru_5/strided_slice_1/stack:output:0&gru_5/strided_slice_1/stack_1:output:0&gru_5/strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
gru_5/strided_slice_1
!gru_5/TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2#
!gru_5/TensorArrayV2/element_shapeЪ
gru_5/TensorArrayV2TensorListReserve*gru_5/TensorArrayV2/element_shape:output:0gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_5/TensorArrayV2Ы
;gru_5/TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2=
;gru_5/TensorArrayUnstack/TensorListFromTensor/element_shape
-gru_5/TensorArrayUnstack/TensorListFromTensorTensorListFromTensorgru_5/transpose:y:0Dgru_5/TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02/
-gru_5/TensorArrayUnstack/TensorListFromTensor
gru_5/strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_2/stack
gru_5/strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_2/stack_1
gru_5/strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_2/stack_2 
gru_5/strided_slice_2StridedSlicegru_5/transpose:y:0$gru_5/strided_slice_2/stack:output:0&gru_5/strided_slice_2/stack_1:output:0&gru_5/strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_5/strided_slice_2
 gru_5/gru_cell_5/ones_like/ShapeShapegru_5/strided_slice_2:output:0*
T0*
_output_shapes
:2"
 gru_5/gru_cell_5/ones_like/Shape
 gru_5/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 gru_5/gru_cell_5/ones_like/ConstШ
gru_5/gru_cell_5/ones_likeFill)gru_5/gru_cell_5/ones_like/Shape:output:0)gru_5/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/ones_like
gru_5/gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
gru_5/gru_cell_5/dropout/ConstУ
gru_5/gru_cell_5/dropout/MulMul#gru_5/gru_cell_5/ones_like:output:0'gru_5/gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/dropout/Mul
gru_5/gru_cell_5/dropout/ShapeShape#gru_5/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2 
gru_5/gru_cell_5/dropout/Shape
5gru_5/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform'gru_5/gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2хЌ27
5gru_5/gru_cell_5/dropout/random_uniform/RandomUniform
'gru_5/gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'gru_5/gru_cell_5/dropout/GreaterEqual/y
%gru_5/gru_cell_5/dropout/GreaterEqualGreaterEqual>gru_5/gru_cell_5/dropout/random_uniform/RandomUniform:output:00gru_5/gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%gru_5/gru_cell_5/dropout/GreaterEqualВ
gru_5/gru_cell_5/dropout/CastCast)gru_5/gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/dropout/CastО
gru_5/gru_cell_5/dropout/Mul_1Mul gru_5/gru_cell_5/dropout/Mul:z:0!gru_5/gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_5/gru_cell_5/dropout/Mul_1
 gru_5/gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 gru_5/gru_cell_5/dropout_1/ConstЩ
gru_5/gru_cell_5/dropout_1/MulMul#gru_5/gru_cell_5/ones_like:output:0)gru_5/gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_5/gru_cell_5/dropout_1/Mul
 gru_5/gru_cell_5/dropout_1/ShapeShape#gru_5/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 gru_5/gru_cell_5/dropout_1/Shape
7gru_5/gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform)gru_5/gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЊЖ29
7gru_5/gru_cell_5/dropout_1/random_uniform/RandomUniform
)gru_5/gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)gru_5/gru_cell_5/dropout_1/GreaterEqual/y
'gru_5/gru_cell_5/dropout_1/GreaterEqualGreaterEqual@gru_5/gru_cell_5/dropout_1/random_uniform/RandomUniform:output:02gru_5/gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'gru_5/gru_cell_5/dropout_1/GreaterEqualИ
gru_5/gru_cell_5/dropout_1/CastCast+gru_5/gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
gru_5/gru_cell_5/dropout_1/CastЦ
 gru_5/gru_cell_5/dropout_1/Mul_1Mul"gru_5/gru_cell_5/dropout_1/Mul:z:0#gru_5/gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/gru_cell_5/dropout_1/Mul_1
 gru_5/gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 gru_5/gru_cell_5/dropout_2/ConstЩ
gru_5/gru_cell_5/dropout_2/MulMul#gru_5/gru_cell_5/ones_like:output:0)gru_5/gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_5/gru_cell_5/dropout_2/Mul
 gru_5/gru_cell_5/dropout_2/ShapeShape#gru_5/gru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2"
 gru_5/gru_cell_5/dropout_2/Shape
7gru_5/gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform)gru_5/gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЖЁ29
7gru_5/gru_cell_5/dropout_2/random_uniform/RandomUniform
)gru_5/gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)gru_5/gru_cell_5/dropout_2/GreaterEqual/y
'gru_5/gru_cell_5/dropout_2/GreaterEqualGreaterEqual@gru_5/gru_cell_5/dropout_2/random_uniform/RandomUniform:output:02gru_5/gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'gru_5/gru_cell_5/dropout_2/GreaterEqualИ
gru_5/gru_cell_5/dropout_2/CastCast+gru_5/gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
gru_5/gru_cell_5/dropout_2/CastЦ
 gru_5/gru_cell_5/dropout_2/Mul_1Mul"gru_5/gru_cell_5/dropout_2/Mul:z:0#gru_5/gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/gru_cell_5/dropout_2/Mul_1Ћ
gru_5/gru_cell_5/ReadVariableOpReadVariableOp(gru_5_gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02!
gru_5/gru_cell_5/ReadVariableOp
gru_5/gru_cell_5/unstackUnpack'gru_5/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_5/gru_cell_5/unstackЉ
gru_5/gru_cell_5/mulMulgru_5/strided_slice_2:output:0"gru_5/gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mulР
&gru_5/gru_cell_5/MatMul/ReadVariableOpReadVariableOp/gru_5_gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02(
&gru_5/gru_cell_5/MatMul/ReadVariableOpИ
gru_5/gru_cell_5/MatMulMatMulgru_5/gru_cell_5/mul:z:0.gru_5/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/MatMulЗ
gru_5/gru_cell_5/BiasAddBiasAdd!gru_5/gru_cell_5/MatMul:product:0!gru_5/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/BiasAddr
gru_5/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/gru_cell_5/Const
 gru_5/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 gru_5/gru_cell_5/split/split_dim№
gru_5/gru_cell_5/splitSplit)gru_5/gru_cell_5/split/split_dim:output:0!gru_5/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/gru_cell_5/splitЦ
(gru_5/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02*
(gru_5/gru_cell_5/MatMul_1/ReadVariableOpК
gru_5/gru_cell_5/MatMul_1MatMulgru_5/zeros:output:00gru_5/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/MatMul_1Н
gru_5/gru_cell_5/BiasAdd_1BiasAdd#gru_5/gru_cell_5/MatMul_1:product:0!gru_5/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/gru_cell_5/BiasAdd_1
gru_5/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_5/gru_cell_5/Const_1
"gru_5/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"gru_5/gru_cell_5/split_1/split_dimЈ
gru_5/gru_cell_5/split_1SplitV#gru_5/gru_cell_5/BiasAdd_1:output:0!gru_5/gru_cell_5/Const_1:output:0+gru_5/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/gru_cell_5/split_1Ћ
gru_5/gru_cell_5/addAddV2gru_5/gru_cell_5/split:output:0!gru_5/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add
gru_5/gru_cell_5/SigmoidSigmoidgru_5/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/SigmoidЏ
gru_5/gru_cell_5/add_1AddV2gru_5/gru_cell_5/split:output:1!gru_5/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_1
gru_5/gru_cell_5/Sigmoid_1Sigmoidgru_5/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/Sigmoid_1Ќ
gru_5/gru_cell_5/mul_1Mulgru_5/gru_cell_5/Sigmoid_1:y:0!gru_5/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_1Ј
gru_5/gru_cell_5/add_2AddV2gru_5/gru_cell_5/split:output:2gru_5/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_2
gru_5/gru_cell_5/Sigmoid_2Sigmoidgru_5/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/Sigmoid_2
gru_5/gru_cell_5/mul_2Mulgru_5/gru_cell_5/Sigmoid:y:0gru_5/zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_2u
gru_5/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_5/gru_cell_5/sub/xЄ
gru_5/gru_cell_5/subSubgru_5/gru_cell_5/sub/x:output:0gru_5/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/subЃ
gru_5/gru_cell_5/mul_3Mulgru_5/gru_cell_5/sub:z:0gru_5/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/mul_3Ѓ
gru_5/gru_cell_5/add_3AddV2gru_5/gru_cell_5/mul_2:z:0gru_5/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/gru_cell_5/add_3
#gru_5/TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2%
#gru_5/TensorArrayV2_1/element_shapeа
gru_5/TensorArrayV2_1TensorListReserve,gru_5/TensorArrayV2_1/element_shape:output:0gru_5/strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
gru_5/TensorArrayV2_1Z

gru_5/timeConst*
_output_shapes
: *
dtype0*
value	B : 2

gru_5/time
gru_5/while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2 
gru_5/while/maximum_iterationsv
gru_5/while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
gru_5/while/loop_counterљ
gru_5/whileWhile!gru_5/while/loop_counter:output:0'gru_5/while/maximum_iterations:output:0gru_5/time:output:0gru_5/TensorArrayV2_1:handle:0gru_5/zeros:output:0gru_5/strided_slice_1:output:0=gru_5/TensorArrayUnstack/TensorListFromTensor:output_handle:0(gru_5_gru_cell_5_readvariableop_resource/gru_5_gru_cell_5_matmul_readvariableop_resource1gru_5_gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*"
bodyR
gru_5_while_body_70564*"
condR
gru_5_while_cond_70563*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
gru_5/whileС
6gru_5/TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   28
6gru_5/TensorArrayV2Stack/TensorListStack/element_shape
(gru_5/TensorArrayV2Stack/TensorListStackTensorListStackgru_5/while:output:3?gru_5/TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02*
(gru_5/TensorArrayV2Stack/TensorListStack
gru_5/strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
gru_5/strided_slice_3/stack
gru_5/strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
gru_5/strided_slice_3/stack_1
gru_5/strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
gru_5/strided_slice_3/stack_2О
gru_5/strided_slice_3StridedSlice1gru_5/TensorArrayV2Stack/TensorListStack:tensor:0$gru_5/strided_slice_3/stack:output:0&gru_5/strided_slice_3/stack_1:output:0&gru_5/strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
gru_5/strided_slice_3
gru_5/transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
gru_5/transpose_1/permО
gru_5/transpose_1	Transpose1gru_5/TensorArrayV2Stack/TensorListStack:tensor:0gru_5/transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
gru_5/transpose_1r
gru_5/runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2
gru_5/runtimeЅ
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes

:*
dtype02
dense_2/MatMul/ReadVariableOpЃ
dense_2/MatMulMatMulgru_5/strided_slice_3:output:0%dense_2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/MatMulЄ
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_2/BiasAdd/ReadVariableOpЁ
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/BiasAddy
dense_2/SigmoidSigmoiddense_2/BiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dense_2/Sigmoid
IdentityIdentitydense_2/Sigmoid:y:0^gru_4/while^gru_5/while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2
gru_4/whilegru_4/while2
gru_5/whilegru_5/while:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ы
Ѕ
while_cond_71786
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_71786___redundant_placeholder03
/while_while_cond_71786___redundant_placeholder13
/while_while_cond_71786___redundant_placeholder23
/while_while_cond_71786___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
А 
щ
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72983

inputs
states_0
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likex
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack_
mulMulinputsones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOp{
MatMul_1MatMulstates_0MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2^
mul_2MulSigmoid:y:0states_0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
Х
ъ
,__inference_sequential_2_layer_call_fn_71078

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
identityЂStatefulPartitionedCallв
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_701962
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Ы
Ѕ
while_cond_72403
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_72403___redundant_placeholder03
/while_while_cond_72403___redundant_placeholder13
/while_while_cond_72403___redundant_placeholder23
/while_while_cond_72403___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
ђD
Ў
while_body_72191
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackМ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 


%__inference_gru_5_layer_call_fn_72307
inputs_0
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCallџ
StatefulPartitionedCallStatefulPartitionedCallinputs_0unknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_692072
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
;
щ
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72939

inputs
states_0
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likec
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout/Const
dropout/MulMulones_like:output:0dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul`
dropout/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout/Shapeг
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Оиу2&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout/GreaterEqual/yО
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/GreaterEqual
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout/Castz
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout/Mul_1g
dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_1/Const
dropout_1/MulMulones_like:output:0dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Muld
dropout_1/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_1/Shapeй
&dropout_1/random_uniform/RandomUniformRandomUniformdropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2дТш2(
&dropout_1/random_uniform/RandomUniformy
dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_1/GreaterEqual/yЦ
dropout_1/GreaterEqualGreaterEqual/dropout_1/random_uniform/RandomUniform:output:0!dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/GreaterEqual
dropout_1/CastCastdropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Cast
dropout_1/Mul_1Muldropout_1/Mul:z:0dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_1/Mul_1g
dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
dropout_2/Const
dropout_2/MulMulones_like:output:0dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Muld
dropout_2/ShapeShapeones_like:output:0*
T0*
_output_shapes
:2
dropout_2/Shapeй
&dropout_2/random_uniform/RandomUniformRandomUniformdropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ќК2(
&dropout_2/random_uniform/RandomUniformy
dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2
dropout_2/GreaterEqual/yЦ
dropout_2/GreaterEqualGreaterEqual/dropout_2/random_uniform/RandomUniform:output:0!dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/GreaterEqual
dropout_2/CastCastdropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Cast
dropout_2/Mul_1Muldropout_2/Mul:z:0dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
dropout_2/Mul_1x
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack^
mulMulinputsdropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOp{
MatMul_1MatMulstates_0MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2^
mul_2MulSigmoid:y:0states_0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
к	
Ќ
*__inference_gru_cell_5_layer_call_fn_72997

inputs
states_0
unknown
	unknown_0
	unknown_1
identity

identity_1ЂStatefulPartitionedCallЁ
StatefulPartitionedCallStatefulPartitionedCallinputsstates_0unknown	unknown_0	unknown_1*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687222
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
Ј 
ч
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_68172

inputs

states
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likex
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack_
mulMulinputsones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOpy
MatMul_1MatMulstatesMatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2\
mul_2MulSigmoid:y:0states*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:OK
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_namestates
Љ
Њ
B__inference_dense_2_layer_call_and_return_conditional_losses_72722

inputs"
matmul_readvariableop_resource#
biasadd_readvariableop_resource
identity
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2
MatMul
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoid_
IdentityIdentitySigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*.
_input_shapes
:џџџџџџџџџ:::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs
фa
Џ

#sequential_2_gru_5_while_body_67927B
>sequential_2_gru_5_while_sequential_2_gru_5_while_loop_counterH
Dsequential_2_gru_5_while_sequential_2_gru_5_while_maximum_iterations(
$sequential_2_gru_5_while_placeholder*
&sequential_2_gru_5_while_placeholder_1*
&sequential_2_gru_5_while_placeholder_2A
=sequential_2_gru_5_while_sequential_2_gru_5_strided_slice_1_0}
ysequential_2_gru_5_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_5_tensorarrayunstack_tensorlistfromtensor_0A
=sequential_2_gru_5_while_gru_cell_5_readvariableop_resource_0H
Dsequential_2_gru_5_while_gru_cell_5_matmul_readvariableop_resource_0J
Fsequential_2_gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0%
!sequential_2_gru_5_while_identity'
#sequential_2_gru_5_while_identity_1'
#sequential_2_gru_5_while_identity_2'
#sequential_2_gru_5_while_identity_3'
#sequential_2_gru_5_while_identity_4?
;sequential_2_gru_5_while_sequential_2_gru_5_strided_slice_1{
wsequential_2_gru_5_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_5_tensorarrayunstack_tensorlistfromtensor?
;sequential_2_gru_5_while_gru_cell_5_readvariableop_resourceF
Bsequential_2_gru_5_while_gru_cell_5_matmul_readvariableop_resourceH
Dsequential_2_gru_5_while_gru_cell_5_matmul_1_readvariableop_resourceщ
Jsequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2L
Jsequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeХ
<sequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItemTensorListGetItemysequential_2_gru_5_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_5_tensorarrayunstack_tensorlistfromtensor_0$sequential_2_gru_5_while_placeholderSsequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02>
<sequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItemн
3sequential_2/gru_5/while/gru_cell_5/ones_like/ShapeShapeCsequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:25
3sequential_2/gru_5/while/gru_cell_5/ones_like/ShapeЏ
3sequential_2/gru_5/while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?25
3sequential_2/gru_5/while/gru_cell_5/ones_like/Const
-sequential_2/gru_5/while/gru_cell_5/ones_likeFill<sequential_2/gru_5/while/gru_cell_5/ones_like/Shape:output:0<sequential_2/gru_5/while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_5/while/gru_cell_5/ones_likeц
2sequential_2/gru_5/while/gru_cell_5/ReadVariableOpReadVariableOp=sequential_2_gru_5_while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype024
2sequential_2/gru_5/while/gru_cell_5/ReadVariableOpж
+sequential_2/gru_5/while/gru_cell_5/unstackUnpack:sequential_2/gru_5/while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2-
+sequential_2/gru_5/while/gru_cell_5/unstack
'sequential_2/gru_5/while/gru_cell_5/mulMulCsequential_2/gru_5/while/TensorArrayV2Read/TensorListGetItem:item:06sequential_2/gru_5/while/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/while/gru_cell_5/mulћ
9sequential_2/gru_5/while/gru_cell_5/MatMul/ReadVariableOpReadVariableOpDsequential_2_gru_5_while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02;
9sequential_2/gru_5/while/gru_cell_5/MatMul/ReadVariableOp
*sequential_2/gru_5/while/gru_cell_5/MatMulMatMul+sequential_2/gru_5/while/gru_cell_5/mul:z:0Asequential_2/gru_5/while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02,
*sequential_2/gru_5/while/gru_cell_5/MatMul
+sequential_2/gru_5/while/gru_cell_5/BiasAddBiasAdd4sequential_2/gru_5/while/gru_cell_5/MatMul:product:04sequential_2/gru_5/while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02-
+sequential_2/gru_5/while/gru_cell_5/BiasAdd
)sequential_2/gru_5/while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2+
)sequential_2/gru_5/while/gru_cell_5/ConstЕ
3sequential_2/gru_5/while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ25
3sequential_2/gru_5/while/gru_cell_5/split/split_dimМ
)sequential_2/gru_5/while/gru_cell_5/splitSplit<sequential_2/gru_5/while/gru_cell_5/split/split_dim:output:04sequential_2/gru_5/while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2+
)sequential_2/gru_5/while/gru_cell_5/split
;sequential_2/gru_5/while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOpFsequential_2_gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02=
;sequential_2/gru_5/while/gru_cell_5/MatMul_1/ReadVariableOp
,sequential_2/gru_5/while/gru_cell_5/MatMul_1MatMul&sequential_2_gru_5_while_placeholder_2Csequential_2/gru_5/while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02.
,sequential_2/gru_5/while/gru_cell_5/MatMul_1
-sequential_2/gru_5/while/gru_cell_5/BiasAdd_1BiasAdd6sequential_2/gru_5/while/gru_cell_5/MatMul_1:product:04sequential_2/gru_5/while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02/
-sequential_2/gru_5/while/gru_cell_5/BiasAdd_1Џ
+sequential_2/gru_5/while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2-
+sequential_2/gru_5/while/gru_cell_5/Const_1Й
5sequential_2/gru_5/while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ27
5sequential_2/gru_5/while/gru_cell_5/split_1/split_dim
+sequential_2/gru_5/while/gru_cell_5/split_1SplitV6sequential_2/gru_5/while/gru_cell_5/BiasAdd_1:output:04sequential_2/gru_5/while/gru_cell_5/Const_1:output:0>sequential_2/gru_5/while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2-
+sequential_2/gru_5/while/gru_cell_5/split_1ї
'sequential_2/gru_5/while/gru_cell_5/addAddV22sequential_2/gru_5/while/gru_cell_5/split:output:04sequential_2/gru_5/while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/while/gru_cell_5/addФ
+sequential_2/gru_5/while/gru_cell_5/SigmoidSigmoid+sequential_2/gru_5/while/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+sequential_2/gru_5/while/gru_cell_5/Sigmoidћ
)sequential_2/gru_5/while/gru_cell_5/add_1AddV22sequential_2/gru_5/while/gru_cell_5/split:output:14sequential_2/gru_5/while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/add_1Ъ
-sequential_2/gru_5/while/gru_cell_5/Sigmoid_1Sigmoid-sequential_2/gru_5/while/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_5/while/gru_cell_5/Sigmoid_1ј
)sequential_2/gru_5/while/gru_cell_5/mul_1Mul1sequential_2/gru_5/while/gru_cell_5/Sigmoid_1:y:04sequential_2/gru_5/while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/mul_1є
)sequential_2/gru_5/while/gru_cell_5/add_2AddV22sequential_2/gru_5/while/gru_cell_5/split:output:2-sequential_2/gru_5/while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/add_2Ъ
-sequential_2/gru_5/while/gru_cell_5/Sigmoid_2Sigmoid-sequential_2/gru_5/while/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_5/while/gru_cell_5/Sigmoid_2ш
)sequential_2/gru_5/while/gru_cell_5/mul_2Mul/sequential_2/gru_5/while/gru_cell_5/Sigmoid:y:0&sequential_2_gru_5_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/mul_2
)sequential_2/gru_5/while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2+
)sequential_2/gru_5/while/gru_cell_5/sub/x№
'sequential_2/gru_5/while/gru_cell_5/subSub2sequential_2/gru_5/while/gru_cell_5/sub/x:output:0/sequential_2/gru_5/while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_5/while/gru_cell_5/subя
)sequential_2/gru_5/while/gru_cell_5/mul_3Mul+sequential_2/gru_5/while/gru_cell_5/sub:z:01sequential_2/gru_5/while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/mul_3я
)sequential_2/gru_5/while/gru_cell_5/add_3AddV2-sequential_2/gru_5/while/gru_cell_5/mul_2:z:0-sequential_2/gru_5/while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_5/while/gru_cell_5/add_3Н
=sequential_2/gru_5/while/TensorArrayV2Write/TensorListSetItemTensorListSetItem&sequential_2_gru_5_while_placeholder_1$sequential_2_gru_5_while_placeholder-sequential_2/gru_5/while/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02?
=sequential_2/gru_5/while/TensorArrayV2Write/TensorListSetItem
sequential_2/gru_5/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2 
sequential_2/gru_5/while/add/yЕ
sequential_2/gru_5/while/addAddV2$sequential_2_gru_5_while_placeholder'sequential_2/gru_5/while/add/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_5/while/add
 sequential_2/gru_5/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2"
 sequential_2/gru_5/while/add_1/yе
sequential_2/gru_5/while/add_1AddV2>sequential_2_gru_5_while_sequential_2_gru_5_while_loop_counter)sequential_2/gru_5/while/add_1/y:output:0*
T0*
_output_shapes
: 2 
sequential_2/gru_5/while/add_1
!sequential_2/gru_5/while/IdentityIdentity"sequential_2/gru_5/while/add_1:z:0*
T0*
_output_shapes
: 2#
!sequential_2/gru_5/while/IdentityН
#sequential_2/gru_5/while/Identity_1IdentityDsequential_2_gru_5_while_sequential_2_gru_5_while_maximum_iterations*
T0*
_output_shapes
: 2%
#sequential_2/gru_5/while/Identity_1
#sequential_2/gru_5/while/Identity_2Identity sequential_2/gru_5/while/add:z:0*
T0*
_output_shapes
: 2%
#sequential_2/gru_5/while/Identity_2Ц
#sequential_2/gru_5/while/Identity_3IdentityMsequential_2/gru_5/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2%
#sequential_2/gru_5/while/Identity_3З
#sequential_2/gru_5/while/Identity_4Identity-sequential_2/gru_5/while/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_5/while/Identity_4"
Dsequential_2_gru_5_while_gru_cell_5_matmul_1_readvariableop_resourceFsequential_2_gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0"
Bsequential_2_gru_5_while_gru_cell_5_matmul_readvariableop_resourceDsequential_2_gru_5_while_gru_cell_5_matmul_readvariableop_resource_0"|
;sequential_2_gru_5_while_gru_cell_5_readvariableop_resource=sequential_2_gru_5_while_gru_cell_5_readvariableop_resource_0"O
!sequential_2_gru_5_while_identity*sequential_2/gru_5/while/Identity:output:0"S
#sequential_2_gru_5_while_identity_1,sequential_2/gru_5/while/Identity_1:output:0"S
#sequential_2_gru_5_while_identity_2,sequential_2/gru_5/while/Identity_2:output:0"S
#sequential_2_gru_5_while_identity_3,sequential_2/gru_5/while/Identity_3:output:0"S
#sequential_2_gru_5_while_identity_4,sequential_2/gru_5/while/Identity_4:output:0"|
;sequential_2_gru_5_while_sequential_2_gru_5_strided_slice_1=sequential_2_gru_5_while_sequential_2_gru_5_strided_slice_1_0"є
wsequential_2_gru_5_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_5_tensorarrayunstack_tensorlistfromtensorysequential_2_gru_5_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_5_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
N
а
gru_5_while_body_70931(
$gru_5_while_gru_5_while_loop_counter.
*gru_5_while_gru_5_while_maximum_iterations
gru_5_while_placeholder
gru_5_while_placeholder_1
gru_5_while_placeholder_2'
#gru_5_while_gru_5_strided_slice_1_0c
_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_04
0gru_5_while_gru_cell_5_readvariableop_resource_0;
7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0=
9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0
gru_5_while_identity
gru_5_while_identity_1
gru_5_while_identity_2
gru_5_while_identity_3
gru_5_while_identity_4%
!gru_5_while_gru_5_strided_slice_1a
]gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor2
.gru_5_while_gru_cell_5_readvariableop_resource9
5gru_5_while_gru_cell_5_matmul_readvariableop_resource;
7gru_5_while_gru_cell_5_matmul_1_readvariableop_resourceЯ
=gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2?
=gru_5/while/TensorArrayV2Read/TensorListGetItem/element_shapeї
/gru_5/while/TensorArrayV2Read/TensorListGetItemTensorListGetItem_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_0gru_5_while_placeholderFgru_5/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype021
/gru_5/while/TensorArrayV2Read/TensorListGetItemЖ
&gru_5/while/gru_cell_5/ones_like/ShapeShape6gru_5/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2(
&gru_5/while/gru_cell_5/ones_like/Shape
&gru_5/while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2(
&gru_5/while/gru_cell_5/ones_like/Constр
 gru_5/while/gru_cell_5/ones_likeFill/gru_5/while/gru_cell_5/ones_like/Shape:output:0/gru_5/while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/ones_likeП
%gru_5/while/gru_cell_5/ReadVariableOpReadVariableOp0gru_5_while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02'
%gru_5/while/gru_cell_5/ReadVariableOpЏ
gru_5/while/gru_cell_5/unstackUnpack-gru_5/while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2 
gru_5/while/gru_cell_5/unstackд
gru_5/while/gru_cell_5/mulMul6gru_5/while/TensorArrayV2Read/TensorListGetItem:item:0)gru_5/while/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mulд
,gru_5/while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02.
,gru_5/while/gru_cell_5/MatMul/ReadVariableOpа
gru_5/while/gru_cell_5/MatMulMatMulgru_5/while/gru_cell_5/mul:z:04gru_5/while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_5/while/gru_cell_5/MatMulЯ
gru_5/while/gru_cell_5/BiasAddBiasAdd'gru_5/while/gru_cell_5/MatMul:product:0'gru_5/while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02 
gru_5/while/gru_cell_5/BiasAdd~
gru_5/while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/gru_cell_5/Const
&gru_5/while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2(
&gru_5/while/gru_cell_5/split/split_dim
gru_5/while/gru_cell_5/splitSplit/gru_5/while/gru_cell_5/split/split_dim:output:0'gru_5/while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_5/while/gru_cell_5/splitк
.gru_5/while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype020
.gru_5/while/gru_cell_5/MatMul_1/ReadVariableOpб
gru_5/while/gru_cell_5/MatMul_1MatMulgru_5_while_placeholder_26gru_5/while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02!
gru_5/while/gru_cell_5/MatMul_1е
 gru_5/while/gru_cell_5/BiasAdd_1BiasAdd)gru_5/while/gru_cell_5/MatMul_1:product:0'gru_5/while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02"
 gru_5/while/gru_cell_5/BiasAdd_1
gru_5/while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2 
gru_5/while/gru_cell_5/Const_1
(gru_5/while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2*
(gru_5/while/gru_cell_5/split_1/split_dimЦ
gru_5/while/gru_cell_5/split_1SplitV)gru_5/while/gru_cell_5/BiasAdd_1:output:0'gru_5/while/gru_cell_5/Const_1:output:01gru_5/while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2 
gru_5/while/gru_cell_5/split_1У
gru_5/while/gru_cell_5/addAddV2%gru_5/while/gru_cell_5/split:output:0'gru_5/while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add
gru_5/while/gru_cell_5/SigmoidSigmoidgru_5/while/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_5/while/gru_cell_5/SigmoidЧ
gru_5/while/gru_cell_5/add_1AddV2%gru_5/while/gru_cell_5/split:output:1'gru_5/while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_1Ѓ
 gru_5/while/gru_cell_5/Sigmoid_1Sigmoid gru_5/while/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/Sigmoid_1Ф
gru_5/while/gru_cell_5/mul_1Mul$gru_5/while/gru_cell_5/Sigmoid_1:y:0'gru_5/while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_1Р
gru_5/while/gru_cell_5/add_2AddV2%gru_5/while/gru_cell_5/split:output:2 gru_5/while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_2Ѓ
 gru_5/while/gru_cell_5/Sigmoid_2Sigmoid gru_5/while/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_5/while/gru_cell_5/Sigmoid_2Д
gru_5/while/gru_cell_5/mul_2Mul"gru_5/while/gru_cell_5/Sigmoid:y:0gru_5_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_2
gru_5/while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_5/while/gru_cell_5/sub/xМ
gru_5/while/gru_cell_5/subSub%gru_5/while/gru_cell_5/sub/x:output:0"gru_5/while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/subЛ
gru_5/while/gru_cell_5/mul_3Mulgru_5/while/gru_cell_5/sub:z:0$gru_5/while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/mul_3Л
gru_5/while/gru_cell_5/add_3AddV2 gru_5/while/gru_cell_5/mul_2:z:0 gru_5/while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/gru_cell_5/add_3ќ
0gru_5/while/TensorArrayV2Write/TensorListSetItemTensorListSetItemgru_5_while_placeholder_1gru_5_while_placeholder gru_5/while/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype022
0gru_5/while/TensorArrayV2Write/TensorListSetItemh
gru_5/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/add/y
gru_5/while/addAddV2gru_5_while_placeholdergru_5/while/add/y:output:0*
T0*
_output_shapes
: 2
gru_5/while/addl
gru_5/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_5/while/add_1/y
gru_5/while/add_1AddV2$gru_5_while_gru_5_while_loop_countergru_5/while/add_1/y:output:0*
T0*
_output_shapes
: 2
gru_5/while/add_1p
gru_5/while/IdentityIdentitygru_5/while/add_1:z:0*
T0*
_output_shapes
: 2
gru_5/while/Identity
gru_5/while/Identity_1Identity*gru_5_while_gru_5_while_maximum_iterations*
T0*
_output_shapes
: 2
gru_5/while/Identity_1r
gru_5/while/Identity_2Identitygru_5/while/add:z:0*
T0*
_output_shapes
: 2
gru_5/while/Identity_2
gru_5/while/Identity_3Identity@gru_5/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
gru_5/while/Identity_3
gru_5/while/Identity_4Identity gru_5/while/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_5/while/Identity_4"H
!gru_5_while_gru_5_strided_slice_1#gru_5_while_gru_5_strided_slice_1_0"t
7gru_5_while_gru_cell_5_matmul_1_readvariableop_resource9gru_5_while_gru_cell_5_matmul_1_readvariableop_resource_0"p
5gru_5_while_gru_cell_5_matmul_readvariableop_resource7gru_5_while_gru_cell_5_matmul_readvariableop_resource_0"b
.gru_5_while_gru_cell_5_readvariableop_resource0gru_5_while_gru_cell_5_readvariableop_resource_0"5
gru_5_while_identitygru_5/while/Identity:output:0"9
gru_5_while_identity_1gru_5/while/Identity_1:output:0"9
gru_5_while_identity_2gru_5/while/Identity_2:output:0"9
gru_5_while_identity_3gru_5/while/Identity_3:output:0"9
gru_5_while_identity_4gru_5/while/Identity_4:output:0"Р
]gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_gru_5_while_tensorarrayv2read_tensorlistgetitem_gru_5_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Ј 
ч
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_68766

inputs

states
readvariableop_resource"
matmul_readvariableop_resource$
 matmul_1_readvariableop_resource
identity

identity_1X
ones_like/ShapeShapeinputs*
T0*
_output_shapes
:2
ones_like/Shapeg
ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
ones_like/Const
	ones_likeFillones_like/Shape:output:0ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	ones_likex
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes

:0*
dtype02
ReadVariableOpj
unstackUnpackReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2	
unstack_
mulMulinputsones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulmul:z:0MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
MatMuls
BiasAddBiasAddMatMul:product:0unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02	
BiasAddP
ConstConst*
_output_shapes
: *
dtype0*
value	B :2
Constm
split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split/split_dimЌ
splitSplitsplit/split_dim:output:0BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
split
MatMul_1/ReadVariableOpReadVariableOp matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02
MatMul_1/ReadVariableOpy
MatMul_1MatMulstatesMatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02

MatMul_1y
	BiasAdd_1BiasAddMatMul_1:product:0unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
	BiasAdd_1g
Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2	
Const_1q
split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
split_1/split_dimг
split_1SplitVBiasAdd_1:output:0Const_1:output:0split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2	
split_1g
addAddV2split:output:0split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
addX
SigmoidSigmoidadd:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2	
Sigmoidk
add_1AddV2split:output:1split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
add_1^
	Sigmoid_1Sigmoid	add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_1h
mul_1MulSigmoid_1:y:0split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_1d
add_2AddV2split:output:2	mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_2^
	Sigmoid_2Sigmoid	add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
	Sigmoid_2\
mul_2MulSigmoid:y:0states*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_2S
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
sub/x`
subSubsub/x:output:0Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
sub_
mul_3Mulsub:z:0Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
mul_3_
add_3AddV2	mul_2:z:0	mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
add_3]
IdentityIdentity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identitya

Identity_1Identity	add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ::::O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:OK
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_namestates
С!
Э
while_body_68431
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0
while_gru_cell_4_68453_0
while_gru_cell_4_68455_0
while_gru_cell_4_68457_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor
while_gru_cell_4_68453
while_gru_cell_4_68455
while_gru_cell_4_68457Ђ(while/gru_cell_4/StatefulPartitionedCallУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЇ
(while/gru_cell_4/StatefulPartitionedCallStatefulPartitionedCall0while/TensorArrayV2Read/TensorListGetItem:item:0while_placeholder_2while_gru_cell_4_68453_0while_gru_cell_4_68455_0while_gru_cell_4_68457_0*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681282*
(while/gru_cell_4/StatefulPartitionedCallѕ
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholder1while/gru_cell_4/StatefulPartitionedCall:output:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1
while/IdentityIdentitywhile/add_1:z:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity
while/Identity_1Identitywhile_while_maximum_iterations)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_1
while/Identity_2Identitywhile/add:z:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_2И
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_3Р
while/Identity_4Identity1while/gru_cell_4/StatefulPartitionedCall:output:1)^while/gru_cell_4/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"2
while_gru_cell_4_68453while_gru_cell_4_68453_0"2
while_gru_cell_4_68455while_gru_cell_4_68455_0"2
while_gru_cell_4_68457while_gru_cell_4_68457_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::2T
(while/gru_cell_4/StatefulPartitionedCall(while/gru_cell_4/StatefulPartitionedCall: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
я[
я
@__inference_gru_5_layer_call_and_return_conditional_losses_72689

inputs&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_like
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_72595*
condR
while_cond_72594*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
N
а
gru_4_while_body_70768(
$gru_4_while_gru_4_while_loop_counter.
*gru_4_while_gru_4_while_maximum_iterations
gru_4_while_placeholder
gru_4_while_placeholder_1
gru_4_while_placeholder_2'
#gru_4_while_gru_4_strided_slice_1_0c
_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_04
0gru_4_while_gru_cell_4_readvariableop_resource_0;
7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0=
9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0
gru_4_while_identity
gru_4_while_identity_1
gru_4_while_identity_2
gru_4_while_identity_3
gru_4_while_identity_4%
!gru_4_while_gru_4_strided_slice_1a
]gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor2
.gru_4_while_gru_cell_4_readvariableop_resource9
5gru_4_while_gru_cell_4_matmul_readvariableop_resource;
7gru_4_while_gru_cell_4_matmul_1_readvariableop_resourceЯ
=gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2?
=gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeї
/gru_4/while/TensorArrayV2Read/TensorListGetItemTensorListGetItem_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_0gru_4_while_placeholderFgru_4/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype021
/gru_4/while/TensorArrayV2Read/TensorListGetItemЖ
&gru_4/while/gru_cell_4/ones_like/ShapeShape6gru_4/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2(
&gru_4/while/gru_cell_4/ones_like/Shape
&gru_4/while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2(
&gru_4/while/gru_cell_4/ones_like/Constр
 gru_4/while/gru_cell_4/ones_likeFill/gru_4/while/gru_cell_4/ones_like/Shape:output:0/gru_4/while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/ones_likeП
%gru_4/while/gru_cell_4/ReadVariableOpReadVariableOp0gru_4_while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02'
%gru_4/while/gru_cell_4/ReadVariableOpЏ
gru_4/while/gru_cell_4/unstackUnpack-gru_4/while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2 
gru_4/while/gru_cell_4/unstackд
gru_4/while/gru_cell_4/mulMul6gru_4/while/TensorArrayV2Read/TensorListGetItem:item:0)gru_4/while/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mulд
,gru_4/while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02.
,gru_4/while/gru_cell_4/MatMul/ReadVariableOpа
gru_4/while/gru_cell_4/MatMulMatMulgru_4/while/gru_cell_4/mul:z:04gru_4/while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_4/while/gru_cell_4/MatMulЯ
gru_4/while/gru_cell_4/BiasAddBiasAdd'gru_4/while/gru_cell_4/MatMul:product:0'gru_4/while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02 
gru_4/while/gru_cell_4/BiasAdd~
gru_4/while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/gru_cell_4/Const
&gru_4/while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2(
&gru_4/while/gru_cell_4/split/split_dim
gru_4/while/gru_cell_4/splitSplit/gru_4/while/gru_cell_4/split/split_dim:output:0'gru_4/while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_4/while/gru_cell_4/splitк
.gru_4/while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype020
.gru_4/while/gru_cell_4/MatMul_1/ReadVariableOpб
gru_4/while/gru_cell_4/MatMul_1MatMulgru_4_while_placeholder_26gru_4/while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02!
gru_4/while/gru_cell_4/MatMul_1е
 gru_4/while/gru_cell_4/BiasAdd_1BiasAdd)gru_4/while/gru_cell_4/MatMul_1:product:0'gru_4/while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02"
 gru_4/while/gru_cell_4/BiasAdd_1
gru_4/while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2 
gru_4/while/gru_cell_4/Const_1
(gru_4/while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2*
(gru_4/while/gru_cell_4/split_1/split_dimЦ
gru_4/while/gru_cell_4/split_1SplitV)gru_4/while/gru_cell_4/BiasAdd_1:output:0'gru_4/while/gru_cell_4/Const_1:output:01gru_4/while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2 
gru_4/while/gru_cell_4/split_1У
gru_4/while/gru_cell_4/addAddV2%gru_4/while/gru_cell_4/split:output:0'gru_4/while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add
gru_4/while/gru_cell_4/SigmoidSigmoidgru_4/while/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
gru_4/while/gru_cell_4/SigmoidЧ
gru_4/while/gru_cell_4/add_1AddV2%gru_4/while/gru_cell_4/split:output:1'gru_4/while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_1Ѓ
 gru_4/while/gru_cell_4/Sigmoid_1Sigmoid gru_4/while/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/Sigmoid_1Ф
gru_4/while/gru_cell_4/mul_1Mul$gru_4/while/gru_cell_4/Sigmoid_1:y:0'gru_4/while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_1Р
gru_4/while/gru_cell_4/add_2AddV2%gru_4/while/gru_cell_4/split:output:2 gru_4/while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_2Ѓ
 gru_4/while/gru_cell_4/Sigmoid_2Sigmoid gru_4/while/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 gru_4/while/gru_cell_4/Sigmoid_2Д
gru_4/while/gru_cell_4/mul_2Mul"gru_4/while/gru_cell_4/Sigmoid:y:0gru_4_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_2
gru_4/while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_4/while/gru_cell_4/sub/xМ
gru_4/while/gru_cell_4/subSub%gru_4/while/gru_cell_4/sub/x:output:0"gru_4/while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/subЛ
gru_4/while/gru_cell_4/mul_3Mulgru_4/while/gru_cell_4/sub:z:0$gru_4/while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/mul_3Л
gru_4/while/gru_cell_4/add_3AddV2 gru_4/while/gru_cell_4/mul_2:z:0 gru_4/while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/gru_cell_4/add_3ќ
0gru_4/while/TensorArrayV2Write/TensorListSetItemTensorListSetItemgru_4_while_placeholder_1gru_4_while_placeholder gru_4/while/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype022
0gru_4/while/TensorArrayV2Write/TensorListSetItemh
gru_4/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/add/y
gru_4/while/addAddV2gru_4_while_placeholdergru_4/while/add/y:output:0*
T0*
_output_shapes
: 2
gru_4/while/addl
gru_4/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
gru_4/while/add_1/y
gru_4/while/add_1AddV2$gru_4_while_gru_4_while_loop_countergru_4/while/add_1/y:output:0*
T0*
_output_shapes
: 2
gru_4/while/add_1p
gru_4/while/IdentityIdentitygru_4/while/add_1:z:0*
T0*
_output_shapes
: 2
gru_4/while/Identity
gru_4/while/Identity_1Identity*gru_4_while_gru_4_while_maximum_iterations*
T0*
_output_shapes
: 2
gru_4/while/Identity_1r
gru_4/while/Identity_2Identitygru_4/while/add:z:0*
T0*
_output_shapes
: 2
gru_4/while/Identity_2
gru_4/while/Identity_3Identity@gru_4/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
gru_4/while/Identity_3
gru_4/while/Identity_4Identity gru_4/while/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_4/while/Identity_4"H
!gru_4_while_gru_4_strided_slice_1#gru_4_while_gru_4_strided_slice_1_0"t
7gru_4_while_gru_cell_4_matmul_1_readvariableop_resource9gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0"p
5gru_4_while_gru_cell_4_matmul_readvariableop_resource7gru_4_while_gru_cell_4_matmul_readvariableop_resource_0"b
.gru_4_while_gru_cell_4_readvariableop_resource0gru_4_while_gru_cell_4_readvariableop_resource_0"5
gru_4_while_identitygru_4/while/Identity:output:0"9
gru_4_while_identity_1gru_4/while/Identity_1:output:0"9
gru_4_while_identity_2gru_4/while/Identity_2:output:0"9
gru_4_while_identity_3gru_4/while/Identity_3:output:0"9
gru_4_while_identity_4gru_4/while/Identity_4:output:0"Р
]gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_gru_4_while_tensorarrayv2read_tensorlistgetitem_gru_4_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
њ

gru_5_while_cond_70563(
$gru_5_while_gru_5_while_loop_counter.
*gru_5_while_gru_5_while_maximum_iterations
gru_5_while_placeholder
gru_5_while_placeholder_1
gru_5_while_placeholder_2*
&gru_5_while_less_gru_5_strided_slice_1?
;gru_5_while_gru_5_while_cond_70563___redundant_placeholder0?
;gru_5_while_gru_5_while_cond_70563___redundant_placeholder1?
;gru_5_while_gru_5_while_cond_70563___redundant_placeholder2?
;gru_5_while_gru_5_while_cond_70563___redundant_placeholder3
gru_5_while_identity

gru_5/while/LessLessgru_5_while_placeholder&gru_5_while_less_gru_5_strided_slice_1*
T0*
_output_shapes
: 2
gru_5/while/Lesso
gru_5/while/IdentityIdentitygru_5/while/Less:z:0*
T0
*
_output_shapes
: 2
gru_5/while/Identity"5
gru_5_while_identitygru_5/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Ы
Ѕ
while_cond_72190
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_72190___redundant_placeholder03
/while_while_cond_72190___redundant_placeholder13
/while_while_cond_72190___redundant_placeholder23
/while_while_cond_72190___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:


%__inference_gru_4_layer_call_fn_71499

inputs
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCall
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_696232
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::22
StatefulPartitionedCallStatefulPartitionedCall:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
П

G__inference_sequential_2_layer_call_and_return_conditional_losses_70092
embedding_2_input
embedding_2_69239
gru_4_69646
gru_4_69648
gru_4_69650
gru_5_70057
gru_5_70059
gru_5_70061
dense_2_70086
dense_2_70088
identityЂdense_2/StatefulPartitionedCallЂ#embedding_2/StatefulPartitionedCallЂgru_4/StatefulPartitionedCallЂgru_5/StatefulPartitionedCall
#embedding_2/StatefulPartitionedCallStatefulPartitionedCallembedding_2_inputembedding_2_69239*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*#
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *O
fJRH
F__inference_embedding_2_layer_call_and_return_conditional_losses_692302%
#embedding_2/StatefulPartitionedCallМ
gru_4/StatefulPartitionedCallStatefulPartitionedCall,embedding_2/StatefulPartitionedCall:output:0gru_4_69646gru_4_69648gru_4_69650*
Tin
2*
Tout
2*
_collective_manager_ids
 *,
_output_shapes
:џџџџџџџџџш*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_4_layer_call_and_return_conditional_losses_694562
gru_4/StatefulPartitionedCallБ
gru_5/StatefulPartitionedCallStatefulPartitionedCall&gru_4/StatefulPartitionedCall:output:0gru_5_70057gru_5_70059gru_5_70061*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_698672
gru_5/StatefulPartitionedCallЌ
dense_2/StatefulPartitionedCallStatefulPartitionedCall&gru_5/StatefulPartitionedCall:output:0dense_2_70086dense_2_70088*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *K
fFRD
B__inference_dense_2_layer_call_and_return_conditional_losses_700752!
dense_2/StatefulPartitionedCall
IdentityIdentity(dense_2/StatefulPartitionedCall:output:0 ^dense_2/StatefulPartitionedCall$^embedding_2/StatefulPartitionedCall^gru_4/StatefulPartitionedCall^gru_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2J
#embedding_2/StatefulPartitionedCall#embedding_2/StatefulPartitionedCall2>
gru_4/StatefulPartitionedCallgru_4/StatefulPartitionedCall2>
gru_5/StatefulPartitionedCallgru_5/StatefulPartitionedCall:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input


%__inference_gru_5_layer_call_fn_72296
inputs_0
unknown
	unknown_0
	unknown_1
identityЂStatefulPartitionedCallџ
StatefulPartitionedCallStatefulPartitionedCallinputs_0unknown	unknown_0	unknown_1*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *I
fDRB
@__inference_gru_5_layer_call_and_return_conditional_losses_690892
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
ђD
Ў
while_body_69529
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackМ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
дP
О
__inference__traced_save_73142
file_prefix5
1savev2_embedding_2_embeddings_read_readvariableop-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop(
$savev2_adam_iter_read_readvariableop	*
&savev2_adam_beta_1_read_readvariableop*
&savev2_adam_beta_2_read_readvariableop)
%savev2_adam_decay_read_readvariableop1
-savev2_adam_learning_rate_read_readvariableop6
2savev2_gru_4_gru_cell_4_kernel_read_readvariableop@
<savev2_gru_4_gru_cell_4_recurrent_kernel_read_readvariableop4
0savev2_gru_4_gru_cell_4_bias_read_readvariableop6
2savev2_gru_5_gru_cell_5_kernel_read_readvariableop@
<savev2_gru_5_gru_cell_5_recurrent_kernel_read_readvariableop4
0savev2_gru_5_gru_cell_5_bias_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop<
8savev2_adam_embedding_2_embeddings_m_read_readvariableop4
0savev2_adam_dense_2_kernel_m_read_readvariableop2
.savev2_adam_dense_2_bias_m_read_readvariableop=
9savev2_adam_gru_4_gru_cell_4_kernel_m_read_readvariableopG
Csavev2_adam_gru_4_gru_cell_4_recurrent_kernel_m_read_readvariableop;
7savev2_adam_gru_4_gru_cell_4_bias_m_read_readvariableop=
9savev2_adam_gru_5_gru_cell_5_kernel_m_read_readvariableopG
Csavev2_adam_gru_5_gru_cell_5_recurrent_kernel_m_read_readvariableop;
7savev2_adam_gru_5_gru_cell_5_bias_m_read_readvariableop<
8savev2_adam_embedding_2_embeddings_v_read_readvariableop4
0savev2_adam_dense_2_kernel_v_read_readvariableop2
.savev2_adam_dense_2_bias_v_read_readvariableop=
9savev2_adam_gru_4_gru_cell_4_kernel_v_read_readvariableopG
Csavev2_adam_gru_4_gru_cell_4_recurrent_kernel_v_read_readvariableop;
7savev2_adam_gru_4_gru_cell_4_bias_v_read_readvariableop=
9savev2_adam_gru_5_gru_cell_5_kernel_v_read_readvariableopG
Csavev2_adam_gru_5_gru_cell_5_recurrent_kernel_v_read_readvariableop;
7savev2_adam_gru_5_gru_cell_5_bias_v_read_readvariableop
savev2_const

identity_1ЂMergeV2Checkpoints
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Const
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*<
value3B1 B+_temp_b3c4ddd4702f420d85360303467f9ed4/part2	
Const_1
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shardІ
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilenameи
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:%*
dtype0*ъ
valueрBн%B:layer_with_weights-0/embeddings/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB)optimizer/iter/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_1/.ATTRIBUTES/VARIABLE_VALUEB+optimizer/beta_2/.ATTRIBUTES/VARIABLE_VALUEB*optimizer/decay/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/learning_rate/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/1/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/2/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/3/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/4/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/5/.ATTRIBUTES/VARIABLE_VALUEB0trainable_variables/6/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEBVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/m/.ATTRIBUTES/VARIABLE_VALUEBVlayer_with_weights-0/embeddings/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBRlayer_with_weights-3/kernel/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBPlayer_with_weights-3/bias/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/1/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/2/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/3/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/4/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/5/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEBLtrainable_variables/6/.OPTIMIZER_SLOT/optimizer/v/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_namesв
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:%*
dtype0*]
valueTBR%B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slicesЁ
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:01savev2_embedding_2_embeddings_read_readvariableop)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop$savev2_adam_iter_read_readvariableop&savev2_adam_beta_1_read_readvariableop&savev2_adam_beta_2_read_readvariableop%savev2_adam_decay_read_readvariableop-savev2_adam_learning_rate_read_readvariableop2savev2_gru_4_gru_cell_4_kernel_read_readvariableop<savev2_gru_4_gru_cell_4_recurrent_kernel_read_readvariableop0savev2_gru_4_gru_cell_4_bias_read_readvariableop2savev2_gru_5_gru_cell_5_kernel_read_readvariableop<savev2_gru_5_gru_cell_5_recurrent_kernel_read_readvariableop0savev2_gru_5_gru_cell_5_bias_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop8savev2_adam_embedding_2_embeddings_m_read_readvariableop0savev2_adam_dense_2_kernel_m_read_readvariableop.savev2_adam_dense_2_bias_m_read_readvariableop9savev2_adam_gru_4_gru_cell_4_kernel_m_read_readvariableopCsavev2_adam_gru_4_gru_cell_4_recurrent_kernel_m_read_readvariableop7savev2_adam_gru_4_gru_cell_4_bias_m_read_readvariableop9savev2_adam_gru_5_gru_cell_5_kernel_m_read_readvariableopCsavev2_adam_gru_5_gru_cell_5_recurrent_kernel_m_read_readvariableop7savev2_adam_gru_5_gru_cell_5_bias_m_read_readvariableop8savev2_adam_embedding_2_embeddings_v_read_readvariableop0savev2_adam_dense_2_kernel_v_read_readvariableop.savev2_adam_dense_2_bias_v_read_readvariableop9savev2_adam_gru_4_gru_cell_4_kernel_v_read_readvariableopCsavev2_adam_gru_4_gru_cell_4_recurrent_kernel_v_read_readvariableop7savev2_adam_gru_4_gru_cell_4_bias_v_read_readvariableop9savev2_adam_gru_5_gru_cell_5_kernel_v_read_readvariableopCsavev2_adam_gru_5_gru_cell_5_recurrent_kernel_v_read_readvariableop7savev2_adam_gru_5_gru_cell_5_bias_v_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *3
dtypes)
'2%	2
SaveV2К
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixesЁ
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identitym

Identity_1IdentityIdentity:output:0^MergeV2Checkpoints*
T0*
_output_shapes
: 2

Identity_1"!

identity_1Identity_1:output:0*А
_input_shapes
: :	::: : : : : :0:0:0:0:0:0: : : : :	:::0:0:0:0:0:0:	:::0:0:0:0:0:0: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:%!

_output_shapes
:	:$ 

_output_shapes

:: 

_output_shapes
::

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :$	 

_output_shapes

:0:$
 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :%!

_output_shapes
:	:$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:$ 

_output_shapes

:0:%!

_output_shapes
:	:$ 

_output_shapes

:: 

_output_shapes
::$ 

_output_shapes

:0:$  

_output_shapes

:0:$! 

_output_shapes

:0:$" 

_output_shapes

:0:$# 

_output_shapes

:0:$$ 

_output_shapes

:0:%

_output_shapes
: 
є

#sequential_2_gru_5_while_cond_67926B
>sequential_2_gru_5_while_sequential_2_gru_5_while_loop_counterH
Dsequential_2_gru_5_while_sequential_2_gru_5_while_maximum_iterations(
$sequential_2_gru_5_while_placeholder*
&sequential_2_gru_5_while_placeholder_1*
&sequential_2_gru_5_while_placeholder_2D
@sequential_2_gru_5_while_less_sequential_2_gru_5_strided_slice_1Y
Usequential_2_gru_5_while_sequential_2_gru_5_while_cond_67926___redundant_placeholder0Y
Usequential_2_gru_5_while_sequential_2_gru_5_while_cond_67926___redundant_placeholder1Y
Usequential_2_gru_5_while_sequential_2_gru_5_while_cond_67926___redundant_placeholder2Y
Usequential_2_gru_5_while_sequential_2_gru_5_while_cond_67926___redundant_placeholder3%
!sequential_2_gru_5_while_identity
Я
sequential_2/gru_5/while/LessLess$sequential_2_gru_5_while_placeholder@sequential_2_gru_5_while_less_sequential_2_gru_5_strided_slice_1*
T0*
_output_shapes
: 2
sequential_2/gru_5/while/Less
!sequential_2/gru_5/while/IdentityIdentity!sequential_2/gru_5/while/Less:z:0*
T0
*
_output_shapes
: 2#
!sequential_2/gru_5/while/Identity"O
!sequential_2_gru_5_while_identity*sequential_2/gru_5/while/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
П}
ё
@__inference_gru_4_layer_call_and_return_conditional_losses_71714
inputs_0&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileF
ShapeShapeinputs_0*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputs_0transpose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_likey
gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout/ConstЋ
gru_cell_4/dropout/MulMulgru_cell_4/ones_like:output:0!gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul
gru_cell_4/dropout/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout/Shapeѓ
/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Жь~21
/gru_cell_4/dropout/random_uniform/RandomUniform
!gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_4/dropout/GreaterEqual/yъ
gru_cell_4/dropout/GreaterEqualGreaterEqual8gru_cell_4/dropout/random_uniform/RandomUniform:output:0*gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_4/dropout/GreaterEqual 
gru_cell_4/dropout/CastCast#gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/CastІ
gru_cell_4/dropout/Mul_1Mulgru_cell_4/dropout/Mul:z:0gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul_1}
gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_1/ConstБ
gru_cell_4/dropout_1/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul
gru_cell_4/dropout_1/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_1/Shapeњ
1gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Чг§23
1gru_cell_4/dropout_1/random_uniform/RandomUniform
#gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_1/GreaterEqual/yђ
!gru_cell_4/dropout_1/GreaterEqualGreaterEqual:gru_cell_4/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_1/GreaterEqualІ
gru_cell_4/dropout_1/CastCast%gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/CastЎ
gru_cell_4/dropout_1/Mul_1Mulgru_cell_4/dropout_1/Mul:z:0gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul_1}
gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_2/ConstБ
gru_cell_4/dropout_2/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul
gru_cell_4/dropout_2/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_2/Shapeњ
1gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2л 23
1gru_cell_4/dropout_2/random_uniform/RandomUniform
#gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_2/GreaterEqual/yђ
!gru_cell_4/dropout_2/GreaterEqualGreaterEqual:gru_cell_4/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_2/GreaterEqualІ
gru_cell_4/dropout_2/CastCast%gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/CastЎ
gru_cell_4/dropout_2/Mul_1Mulgru_cell_4/dropout_2/Mul:z:0gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul_1
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_71596*
condR
while_cond_71595*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimex
IdentityIdentitytranspose_1:y:0^while*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2
whilewhile:^ Z
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
"
_user_specified_name
inputs/0
}
я
@__inference_gru_4_layer_call_and_return_conditional_losses_71310

inputs&
"gru_cell_4_readvariableop_resource-
)gru_cell_4_matmul_readvariableop_resource/
+gru_cell_4_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_4/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_4/ones_like/Shape}
gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/ones_like/ConstА
gru_cell_4/ones_likeFill#gru_cell_4/ones_like/Shape:output:0#gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/ones_likey
gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout/ConstЋ
gru_cell_4/dropout/MulMulgru_cell_4/ones_like:output:0!gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul
gru_cell_4/dropout/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout/Shapeѓ
/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ы&21
/gru_cell_4/dropout/random_uniform/RandomUniform
!gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_4/dropout/GreaterEqual/yъ
gru_cell_4/dropout/GreaterEqualGreaterEqual8gru_cell_4/dropout/random_uniform/RandomUniform:output:0*gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_4/dropout/GreaterEqual 
gru_cell_4/dropout/CastCast#gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/CastІ
gru_cell_4/dropout/Mul_1Mulgru_cell_4/dropout/Mul:z:0gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout/Mul_1}
gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_1/ConstБ
gru_cell_4/dropout_1/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul
gru_cell_4/dropout_1/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_1/Shapeњ
1gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ьіб23
1gru_cell_4/dropout_1/random_uniform/RandomUniform
#gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_1/GreaterEqual/yђ
!gru_cell_4/dropout_1/GreaterEqualGreaterEqual:gru_cell_4/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_1/GreaterEqualІ
gru_cell_4/dropout_1/CastCast%gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/CastЎ
gru_cell_4/dropout_1/Mul_1Mulgru_cell_4/dropout_1/Mul:z:0gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_1/Mul_1}
gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_4/dropout_2/ConstБ
gru_cell_4/dropout_2/MulMulgru_cell_4/ones_like:output:0#gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul
gru_cell_4/dropout_2/ShapeShapegru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_4/dropout_2/Shapeњ
1gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ќ23
1gru_cell_4/dropout_2/random_uniform/RandomUniform
#gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_4/dropout_2/GreaterEqual/yђ
!gru_cell_4/dropout_2/GreaterEqualGreaterEqual:gru_cell_4/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_4/dropout_2/GreaterEqualІ
gru_cell_4/dropout_2/CastCast%gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/CastЎ
gru_cell_4/dropout_2/Mul_1Mulgru_cell_4/dropout_2/Mul:z:0gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/dropout_2/Mul_1
gru_cell_4/ReadVariableOpReadVariableOp"gru_cell_4_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_4/ReadVariableOp
gru_cell_4/unstackUnpack!gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_4/unstack
gru_cell_4/mulMulstrided_slice_2:output:0gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mulЎ
 gru_cell_4/MatMul/ReadVariableOpReadVariableOp)gru_cell_4_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_4/MatMul/ReadVariableOp 
gru_cell_4/MatMulMatMulgru_cell_4/mul:z:0(gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul
gru_cell_4/BiasAddBiasAddgru_cell_4/MatMul:product:0gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAddf
gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_4/Const
gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split/split_dimи
gru_cell_4/splitSplit#gru_cell_4/split/split_dim:output:0gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/splitД
"gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_4_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_4/MatMul_1/ReadVariableOpЂ
gru_cell_4/MatMul_1MatMulzeros:output:0*gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/MatMul_1Ѕ
gru_cell_4/BiasAdd_1BiasAddgru_cell_4/MatMul_1:product:0gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_4/BiasAdd_1}
gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_4/Const_1
gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_4/split_1/split_dim
gru_cell_4/split_1SplitVgru_cell_4/BiasAdd_1:output:0gru_cell_4/Const_1:output:0%gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_4/split_1
gru_cell_4/addAddV2gru_cell_4/split:output:0gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/addy
gru_cell_4/SigmoidSigmoidgru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid
gru_cell_4/add_1AddV2gru_cell_4/split:output:1gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_1
gru_cell_4/Sigmoid_1Sigmoidgru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_1
gru_cell_4/mul_1Mulgru_cell_4/Sigmoid_1:y:0gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_1
gru_cell_4/add_2AddV2gru_cell_4/split:output:2gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_2
gru_cell_4/Sigmoid_2Sigmoidgru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/Sigmoid_2
gru_cell_4/mul_2Mulgru_cell_4/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_2i
gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_4/sub/x
gru_cell_4/subSubgru_cell_4/sub/x:output:0gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/sub
gru_cell_4/mul_3Mulgru_cell_4/sub:z:0gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/mul_3
gru_cell_4/add_3AddV2gru_cell_4/mul_2:z:0gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_4/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_4_readvariableop_resource)gru_cell_4_matmul_readvariableop_resource+gru_cell_4_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_71192*
condR
while_cond_71191*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimep
IdentityIdentitytranspose_1:y:0^while*
T0*,
_output_shapes
:џџџџџџџџџш2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
Х
ъ
,__inference_sequential_2_layer_call_fn_71055

inputs
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
identityЂStatefulPartitionedCallв
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_701472
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
ђD
Ў
while_body_72595
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_5_readvariableop_resource_05
1while_gru_cell_5_matmul_readvariableop_resource_07
3while_gru_cell_5_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_5_readvariableop_resource3
/while_gru_cell_5_matmul_readvariableop_resource5
1while_gru_cell_5_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_5/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_5/ones_like/Shape
 while/gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_5/ones_like/ConstШ
while/gru_cell_5/ones_likeFill)while/gru_cell_5/ones_like/Shape:output:0)while/gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/ones_like­
while/gru_cell_5/ReadVariableOpReadVariableOp*while_gru_cell_5_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_5/ReadVariableOp
while/gru_cell_5/unstackUnpack'while/gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_5/unstackМ
while/gru_cell_5/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mulТ
&while/gru_cell_5/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_5_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_5/MatMul/ReadVariableOpИ
while/gru_cell_5/MatMulMatMulwhile/gru_cell_5/mul:z:0.while/gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMulЗ
while/gru_cell_5/BiasAddBiasAdd!while/gru_cell_5/MatMul:product:0!while/gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAddr
while/gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_5/Const
 while/gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_5/split/split_dim№
while/gru_cell_5/splitSplit)while/gru_cell_5/split/split_dim:output:0!while/gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/splitШ
(while/gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_5_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_5/MatMul_1/ReadVariableOpЙ
while/gru_cell_5/MatMul_1MatMulwhile_placeholder_20while/gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/MatMul_1Н
while/gru_cell_5/BiasAdd_1BiasAdd#while/gru_cell_5/MatMul_1:product:0!while/gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_5/BiasAdd_1
while/gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_5/Const_1
"while/gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_5/split_1/split_dimЈ
while/gru_cell_5/split_1SplitV#while/gru_cell_5/BiasAdd_1:output:0!while/gru_cell_5/Const_1:output:0+while/gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_5/split_1Ћ
while/gru_cell_5/addAddV2while/gru_cell_5/split:output:0!while/gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add
while/gru_cell_5/SigmoidSigmoidwhile/gru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/SigmoidЏ
while/gru_cell_5/add_1AddV2while/gru_cell_5/split:output:1!while/gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_1
while/gru_cell_5/Sigmoid_1Sigmoidwhile/gru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_1Ќ
while/gru_cell_5/mul_1Mulwhile/gru_cell_5/Sigmoid_1:y:0!while/gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_1Ј
while/gru_cell_5/add_2AddV2while/gru_cell_5/split:output:2while/gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_2
while/gru_cell_5/Sigmoid_2Sigmoidwhile/gru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/Sigmoid_2
while/gru_cell_5/mul_2Mulwhile/gru_cell_5/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_2u
while/gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_5/sub/xЄ
while/gru_cell_5/subSubwhile/gru_cell_5/sub/x:output:0while/gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/subЃ
while/gru_cell_5/mul_3Mulwhile/gru_cell_5/sub:z:0while/gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/mul_3Ѓ
while/gru_cell_5/add_3AddV2while/gru_cell_5/mul_2:z:0while/gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_5/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_5/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_5/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_5_matmul_1_readvariableop_resource3while_gru_cell_5_matmul_1_readvariableop_resource_0"d
/while_gru_cell_5_matmul_readvariableop_resource1while_gru_cell_5_matmul_readvariableop_resource_0"V
(while_gru_cell_5_readvariableop_resource*while_gru_cell_5_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Ш<
Ю
@__inference_gru_5_layer_call_and_return_conditional_losses_69089

inputs
gru_cell_5_69013
gru_cell_5_69015
gru_cell_5_69017
identityЂ"gru_cell_5/StatefulPartitionedCallЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputstranspose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2ц
"gru_cell_5/StatefulPartitionedCallStatefulPartitionedCallstrided_slice_2:output:0zeros:output:0gru_cell_5_69013gru_cell_5_69015gru_cell_5_69017*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687222$
"gru_cell_5/StatefulPartitionedCall
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterп
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0gru_cell_5_69013gru_cell_5_69015gru_cell_5_69017*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69025*
condR
while_cond_69024*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtime
IdentityIdentitystrided_slice_3:output:0#^gru_cell_5/StatefulPartitionedCall^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2H
"gru_cell_5/StatefulPartitionedCall"gru_cell_5/StatefulPartitionedCall2
whilewhile:\ X
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs
ц
ѕ
,__inference_sequential_2_layer_call_fn_70168
embedding_2_input
unknown
	unknown_0
	unknown_1
	unknown_2
	unknown_3
	unknown_4
	unknown_5
	unknown_6
	unknown_7
identityЂStatefulPartitionedCallн
StatefulPartitionedCallStatefulPartitionedCallembedding_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7*
Tin
2
*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:џџџџџџџџџ*+
_read_only_resource_inputs
		*-
config_proto

CPU

GPU 2J 8 *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_701472
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*K
_input_shapes:
8:џџџџџџџџџш:::::::::22
StatefulPartitionedCallStatefulPartitionedCall:[ W
(
_output_shapes
:џџџџџџџџџш
+
_user_specified_nameembedding_2_input
к	
Ќ
*__inference_gru_cell_5_layer_call_fn_73011

inputs
states_0
unknown
	unknown_0
	unknown_1
identity

identity_1ЂStatefulPartitionedCallЁ
StatefulPartitionedCallStatefulPartitionedCallinputsstates_0unknown	unknown_0	unknown_1*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687662
StatefulPartitionedCall
IdentityIdentity StatefulPartitionedCall:output:0^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity

Identity_1Identity StatefulPartitionedCall:output:1^StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity_1"
identityIdentity:output:0"!

identity_1Identity_1:output:0*E
_input_shapes4
2:џџџџџџџџџ:џџџџџџџџџ:::22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:џџџџџџџџџ
 
_user_specified_nameinputs:QM
'
_output_shapes
:џџџџџџџџџ
"
_user_specified_name
states/0
}
я
@__inference_gru_5_layer_call_and_return_conditional_losses_69867

inputs&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_likey
gru_cell_5/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout/ConstЋ
gru_cell_5/dropout/MulMulgru_cell_5/ones_like:output:0!gru_cell_5/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul
gru_cell_5/dropout/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout/Shapeє
/gru_cell_5/dropout/random_uniform/RandomUniformRandomUniform!gru_cell_5/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Іѓя21
/gru_cell_5/dropout/random_uniform/RandomUniform
!gru_cell_5/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2#
!gru_cell_5/dropout/GreaterEqual/yъ
gru_cell_5/dropout/GreaterEqualGreaterEqual8gru_cell_5/dropout/random_uniform/RandomUniform:output:0*gru_cell_5/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2!
gru_cell_5/dropout/GreaterEqual 
gru_cell_5/dropout/CastCast#gru_cell_5/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/CastІ
gru_cell_5/dropout/Mul_1Mulgru_cell_5/dropout/Mul:z:0gru_cell_5/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout/Mul_1}
gru_cell_5/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_1/ConstБ
gru_cell_5/dropout_1/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul
gru_cell_5/dropout_1/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_1/Shapeњ
1gru_cell_5/dropout_1/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ХЩ23
1gru_cell_5/dropout_1/random_uniform/RandomUniform
#gru_cell_5/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_1/GreaterEqual/yђ
!gru_cell_5/dropout_1/GreaterEqualGreaterEqual:gru_cell_5/dropout_1/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_1/GreaterEqualІ
gru_cell_5/dropout_1/CastCast%gru_cell_5/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/CastЎ
gru_cell_5/dropout_1/Mul_1Mulgru_cell_5/dropout_1/Mul:z:0gru_cell_5/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_1/Mul_1}
gru_cell_5/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2
gru_cell_5/dropout_2/ConstБ
gru_cell_5/dropout_2/MulMulgru_cell_5/ones_like:output:0#gru_cell_5/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul
gru_cell_5/dropout_2/ShapeShapegru_cell_5/ones_like:output:0*
T0*
_output_shapes
:2
gru_cell_5/dropout_2/Shapeљ
1gru_cell_5/dropout_2/random_uniform/RandomUniformRandomUniform#gru_cell_5/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЯЌ
23
1gru_cell_5/dropout_2/random_uniform/RandomUniform
#gru_cell_5/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2%
#gru_cell_5/dropout_2/GreaterEqual/yђ
!gru_cell_5/dropout_2/GreaterEqualGreaterEqual:gru_cell_5/dropout_2/random_uniform/RandomUniform:output:0,gru_cell_5/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2#
!gru_cell_5/dropout_2/GreaterEqualІ
gru_cell_5/dropout_2/CastCast%gru_cell_5/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/CastЎ
gru_cell_5/dropout_2/Mul_1Mulgru_cell_5/dropout_2/Mul:z:0gru_cell_5/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/dropout_2/Mul_1
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69749*
condR
while_cond_69748*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
С!
Э
while_body_69025
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0
while_gru_cell_5_69047_0
while_gru_cell_5_69049_0
while_gru_cell_5_69051_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor
while_gru_cell_5_69047
while_gru_cell_5_69049
while_gru_cell_5_69051Ђ(while/gru_cell_5/StatefulPartitionedCallУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЇ
(while/gru_cell_5/StatefulPartitionedCallStatefulPartitionedCall0while/TensorArrayV2Read/TensorListGetItem:item:0while_placeholder_2while_gru_cell_5_69047_0while_gru_cell_5_69049_0while_gru_cell_5_69051_0*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687222*
(while/gru_cell_5/StatefulPartitionedCallѕ
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholder1while/gru_cell_5/StatefulPartitionedCall:output:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1
while/IdentityIdentitywhile/add_1:z:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity
while/Identity_1Identitywhile_while_maximum_iterations)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_1
while/Identity_2Identitywhile/add:z:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_2И
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_3Р
while/Identity_4Identity1while/gru_cell_5/StatefulPartitionedCall:output:1)^while/gru_cell_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"2
while_gru_cell_5_69047while_gru_cell_5_69047_0"2
while_gru_cell_5_69049while_gru_cell_5_69049_0"2
while_gru_cell_5_69051while_gru_cell_5_69051_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::2T
(while/gru_cell_5/StatefulPartitionedCall(while/gru_cell_5/StatefulPartitionedCall: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
я[
я
@__inference_gru_5_layer_call_and_return_conditional_losses_70034

inputs&
"gru_cell_5_readvariableop_resource-
)gru_cell_5_matmul_readvariableop_resource/
+gru_cell_5_matmul_1_readvariableop_resource
identityЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm{
	transpose	Transposeinputstranspose/perm:output:0*
T0*,
_output_shapes
:шџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2
gru_cell_5/ones_like/ShapeShapestrided_slice_2:output:0*
T0*
_output_shapes
:2
gru_cell_5/ones_like/Shape}
gru_cell_5/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/ones_like/ConstА
gru_cell_5/ones_likeFill#gru_cell_5/ones_like/Shape:output:0#gru_cell_5/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/ones_like
gru_cell_5/ReadVariableOpReadVariableOp"gru_cell_5_readvariableop_resource*
_output_shapes

:0*
dtype02
gru_cell_5/ReadVariableOp
gru_cell_5/unstackUnpack!gru_cell_5/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
gru_cell_5/unstack
gru_cell_5/mulMulstrided_slice_2:output:0gru_cell_5/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mulЎ
 gru_cell_5/MatMul/ReadVariableOpReadVariableOp)gru_cell_5_matmul_readvariableop_resource*
_output_shapes

:0*
dtype02"
 gru_cell_5/MatMul/ReadVariableOp 
gru_cell_5/MatMulMatMulgru_cell_5/mul:z:0(gru_cell_5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul
gru_cell_5/BiasAddBiasAddgru_cell_5/MatMul:product:0gru_cell_5/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAddf
gru_cell_5/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
gru_cell_5/Const
gru_cell_5/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split/split_dimи
gru_cell_5/splitSplit#gru_cell_5/split/split_dim:output:0gru_cell_5/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/splitД
"gru_cell_5/MatMul_1/ReadVariableOpReadVariableOp+gru_cell_5_matmul_1_readvariableop_resource*
_output_shapes

:0*
dtype02$
"gru_cell_5/MatMul_1/ReadVariableOpЂ
gru_cell_5/MatMul_1MatMulzeros:output:0*gru_cell_5/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/MatMul_1Ѕ
gru_cell_5/BiasAdd_1BiasAddgru_cell_5/MatMul_1:product:0gru_cell_5/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
gru_cell_5/BiasAdd_1}
gru_cell_5/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
gru_cell_5/Const_1
gru_cell_5/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
gru_cell_5/split_1/split_dim
gru_cell_5/split_1SplitVgru_cell_5/BiasAdd_1:output:0gru_cell_5/Const_1:output:0%gru_cell_5/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
gru_cell_5/split_1
gru_cell_5/addAddV2gru_cell_5/split:output:0gru_cell_5/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/addy
gru_cell_5/SigmoidSigmoidgru_cell_5/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid
gru_cell_5/add_1AddV2gru_cell_5/split:output:1gru_cell_5/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_1
gru_cell_5/Sigmoid_1Sigmoidgru_cell_5/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_1
gru_cell_5/mul_1Mulgru_cell_5/Sigmoid_1:y:0gru_cell_5/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_1
gru_cell_5/add_2AddV2gru_cell_5/split:output:2gru_cell_5/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_2
gru_cell_5/Sigmoid_2Sigmoidgru_cell_5/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/Sigmoid_2
gru_cell_5/mul_2Mulgru_cell_5/Sigmoid:y:0zeros:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_2i
gru_cell_5/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
gru_cell_5/sub/x
gru_cell_5/subSubgru_cell_5/sub/x:output:0gru_cell_5/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/sub
gru_cell_5/mul_3Mulgru_cell_5/sub:z:0gru_cell_5/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/mul_3
gru_cell_5/add_3AddV2gru_cell_5/mul_2:z:0gru_cell_5/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
gru_cell_5/add_3
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterЅ
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0"gru_cell_5_readvariableop_resource)gru_cell_5_matmul_readvariableop_resource+gru_cell_5_matmul_1_readvariableop_resource*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_69940*
condR
while_cond_69939*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeщ
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*,
_output_shapes
:шџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permІ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*,
_output_shapes
:џџџџџџџџџш2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtimet
IdentityIdentitystrided_slice_3:output:0^while*
T0*'
_output_shapes
:џџџџџџџџџ2

Identity"
identityIdentity:output:0*7
_input_shapes&
$:џџџџџџџџџш:::2
whilewhile:T P
,
_output_shapes
:џџџџџџџџџш
 
_user_specified_nameinputs
С!
Э
while_body_69143
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0
while_gru_cell_5_69165_0
while_gru_cell_5_69167_0
while_gru_cell_5_69169_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor
while_gru_cell_5_69165
while_gru_cell_5_69167
while_gru_cell_5_69169Ђ(while/gru_cell_5/StatefulPartitionedCallУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЇ
(while/gru_cell_5/StatefulPartitionedCallStatefulPartitionedCall0while/TensorArrayV2Read/TensorListGetItem:item:0while_placeholder_2while_gru_cell_5_69165_0while_gru_cell_5_69167_0while_gru_cell_5_69169_0*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_687662*
(while/gru_cell_5/StatefulPartitionedCallѕ
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholder1while/gru_cell_5/StatefulPartitionedCall:output:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1
while/IdentityIdentitywhile/add_1:z:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity
while/Identity_1Identitywhile_while_maximum_iterations)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_1
while/Identity_2Identitywhile/add:z:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_2И
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0)^while/gru_cell_5/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_3Р
while/Identity_4Identity1while/gru_cell_5/StatefulPartitionedCall:output:1)^while/gru_cell_5/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"2
while_gru_cell_5_69165while_gru_cell_5_69165_0"2
while_gru_cell_5_69167while_gru_cell_5_69167_0"2
while_gru_cell_5_69169while_gru_cell_5_69169_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::2T
(while/gru_cell_5/StatefulPartitionedCall(while/gru_cell_5/StatefulPartitionedCall: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Ы
Ѕ
while_cond_71999
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_71999___redundant_placeholder03
/while_while_cond_71999___redundant_placeholder13
/while_while_cond_71999___redundant_placeholder23
/while_while_cond_71999___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
Ы
Ѕ
while_cond_69528
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_less_strided_slice_13
/while_while_cond_69528___redundant_placeholder03
/while_while_cond_69528___redundant_placeholder13
/while_while_cond_69528___redundant_placeholder23
/while_while_cond_69528___redundant_placeholder3
while_identity
p

while/LessLesswhile_placeholderwhile_less_strided_slice_1*
T0*
_output_shapes
: 2

while/Less]
while/IdentityIdentitywhile/Less:z:0*
T0
*
_output_shapes
: 2
while/Identity")
while_identitywhile/Identity:output:0*@
_input_shapes/
-: : : : :џџџџџџџџџ: ::::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
:
ђD
Ў
while_body_71787
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackМ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0#while/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
С!
Э
while_body_68549
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0
while_gru_cell_4_68571_0
while_gru_cell_4_68573_0
while_gru_cell_4_68575_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor
while_gru_cell_4_68571
while_gru_cell_4_68573
while_gru_cell_4_68575Ђ(while/gru_cell_4/StatefulPartitionedCallУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЇ
(while/gru_cell_4/StatefulPartitionedCallStatefulPartitionedCall0while/TensorArrayV2Read/TensorListGetItem:item:0while_placeholder_2while_gru_cell_4_68571_0while_gru_cell_4_68573_0while_gru_cell_4_68575_0*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681722*
(while/gru_cell_4/StatefulPartitionedCallѕ
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholder1while/gru_cell_4/StatefulPartitionedCall:output:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1
while/IdentityIdentitywhile/add_1:z:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity
while/Identity_1Identitywhile_while_maximum_iterations)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_1
while/Identity_2Identitywhile/add:z:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_2И
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0)^while/gru_cell_4/StatefulPartitionedCall*
T0*
_output_shapes
: 2
while/Identity_3Р
while/Identity_4Identity1while/gru_cell_4/StatefulPartitionedCall:output:1)^while/gru_cell_4/StatefulPartitionedCall*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"2
while_gru_cell_4_68571while_gru_cell_4_68571_0"2
while_gru_cell_4_68573while_gru_cell_4_68573_0"2
while_gru_cell_4_68575while_gru_cell_4_68575_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::2T
(while/gru_cell_4/StatefulPartitionedCall(while/gru_cell_4/StatefulPartitionedCall: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
гi
Ў
while_body_71192
while_while_loop_counter"
while_while_maximum_iterations
while_placeholder
while_placeholder_1
while_placeholder_2
while_strided_slice_1_0W
Swhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0.
*while_gru_cell_4_readvariableop_resource_05
1while_gru_cell_4_matmul_readvariableop_resource_07
3while_gru_cell_4_matmul_1_readvariableop_resource_0
while_identity
while_identity_1
while_identity_2
while_identity_3
while_identity_4
while_strided_slice_1U
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor,
(while_gru_cell_4_readvariableop_resource3
/while_gru_cell_4_matmul_readvariableop_resource5
1while_gru_cell_4_matmul_1_readvariableop_resourceУ
7while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   29
7while/TensorArrayV2Read/TensorListGetItem/element_shapeг
)while/TensorArrayV2Read/TensorListGetItemTensorListGetItemSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0while_placeholder@while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02+
)while/TensorArrayV2Read/TensorListGetItemЄ
 while/gru_cell_4/ones_like/ShapeShape0while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/ones_like/Shape
 while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2"
 while/gru_cell_4/ones_like/ConstШ
while/gru_cell_4/ones_likeFill)while/gru_cell_4/ones_like/Shape:output:0)while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/ones_like
while/gru_cell_4/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2 
while/gru_cell_4/dropout/ConstУ
while/gru_cell_4/dropout/MulMul#while/gru_cell_4/ones_like:output:0'while/gru_cell_4/dropout/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/Mul
while/gru_cell_4/dropout/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2 
while/gru_cell_4/dropout/Shape
5while/gru_cell_4/dropout/random_uniform/RandomUniformRandomUniform'while/gru_cell_4/dropout/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЛВ27
5while/gru_cell_4/dropout/random_uniform/RandomUniform
'while/gru_cell_4/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2)
'while/gru_cell_4/dropout/GreaterEqual/y
%while/gru_cell_4/dropout/GreaterEqualGreaterEqual>while/gru_cell_4/dropout/random_uniform/RandomUniform:output:00while/gru_cell_4/dropout/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2'
%while/gru_cell_4/dropout/GreaterEqualВ
while/gru_cell_4/dropout/CastCast)while/gru_cell_4/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/dropout/CastО
while/gru_cell_4/dropout/Mul_1Mul while/gru_cell_4/dropout/Mul:z:0!while/gru_cell_4/dropout/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout/Mul_1
 while/gru_cell_4/dropout_1/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_1/ConstЩ
while/gru_cell_4/dropout_1/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_1/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_1/Mul
 while/gru_cell_4/dropout_1/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_1/Shape
7while/gru_cell_4/dropout_1/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_1/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2Пщ29
7while/gru_cell_4/dropout_1/random_uniform/RandomUniform
)while/gru_cell_4/dropout_1/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_1/GreaterEqual/y
'while/gru_cell_4/dropout_1/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_1/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_1/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_1/GreaterEqualИ
while/gru_cell_4/dropout_1/CastCast+while/gru_cell_4/dropout_1/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_1/CastЦ
 while/gru_cell_4/dropout_1/Mul_1Mul"while/gru_cell_4/dropout_1/Mul:z:0#while/gru_cell_4/dropout_1/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_1/Mul_1
 while/gru_cell_4/dropout_2/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *   @2"
 while/gru_cell_4/dropout_2/ConstЩ
while/gru_cell_4/dropout_2/MulMul#while/gru_cell_4/ones_like:output:0)while/gru_cell_4/dropout_2/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2 
while/gru_cell_4/dropout_2/Mul
 while/gru_cell_4/dropout_2/ShapeShape#while/gru_cell_4/ones_like:output:0*
T0*
_output_shapes
:2"
 while/gru_cell_4/dropout_2/Shape
7while/gru_cell_4/dropout_2/random_uniform/RandomUniformRandomUniform)while/gru_cell_4/dropout_2/Shape:output:0*
T0*'
_output_shapes
:џџџџџџџџџ*
dtype0*
seedБџх)*
seed2ЗЊ29
7while/gru_cell_4/dropout_2/random_uniform/RandomUniform
)while/gru_cell_4/dropout_2/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *   ?2+
)while/gru_cell_4/dropout_2/GreaterEqual/y
'while/gru_cell_4/dropout_2/GreaterEqualGreaterEqual@while/gru_cell_4/dropout_2/random_uniform/RandomUniform:output:02while/gru_cell_4/dropout_2/GreaterEqual/y:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'while/gru_cell_4/dropout_2/GreaterEqualИ
while/gru_cell_4/dropout_2/CastCast+while/gru_cell_4/dropout_2/GreaterEqual:z:0*

DstT0*

SrcT0
*'
_output_shapes
:џџџџџџџџџ2!
while/gru_cell_4/dropout_2/CastЦ
 while/gru_cell_4/dropout_2/Mul_1Mul"while/gru_cell_4/dropout_2/Mul:z:0#while/gru_cell_4/dropout_2/Cast:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2"
 while/gru_cell_4/dropout_2/Mul_1­
while/gru_cell_4/ReadVariableOpReadVariableOp*while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype02!
while/gru_cell_4/ReadVariableOp
while/gru_cell_4/unstackUnpack'while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2
while/gru_cell_4/unstackЛ
while/gru_cell_4/mulMul0while/TensorArrayV2Read/TensorListGetItem:item:0"while/gru_cell_4/dropout/Mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mulТ
&while/gru_cell_4/MatMul/ReadVariableOpReadVariableOp1while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02(
&while/gru_cell_4/MatMul/ReadVariableOpИ
while/gru_cell_4/MatMulMatMulwhile/gru_cell_4/mul:z:0.while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMulЗ
while/gru_cell_4/BiasAddBiasAdd!while/gru_cell_4/MatMul:product:0!while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAddr
while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2
while/gru_cell_4/Const
 while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2"
 while/gru_cell_4/split/split_dim№
while/gru_cell_4/splitSplit)while/gru_cell_4/split/split_dim:output:0!while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/splitШ
(while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOp3while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02*
(while/gru_cell_4/MatMul_1/ReadVariableOpЙ
while/gru_cell_4/MatMul_1MatMulwhile_placeholder_20while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/MatMul_1Н
while/gru_cell_4/BiasAdd_1BiasAdd#while/gru_cell_4/MatMul_1:product:0!while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02
while/gru_cell_4/BiasAdd_1
while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2
while/gru_cell_4/Const_1
"while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2$
"while/gru_cell_4/split_1/split_dimЈ
while/gru_cell_4/split_1SplitV#while/gru_cell_4/BiasAdd_1:output:0!while/gru_cell_4/Const_1:output:0+while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2
while/gru_cell_4/split_1Ћ
while/gru_cell_4/addAddV2while/gru_cell_4/split:output:0!while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add
while/gru_cell_4/SigmoidSigmoidwhile/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/SigmoidЏ
while/gru_cell_4/add_1AddV2while/gru_cell_4/split:output:1!while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_1
while/gru_cell_4/Sigmoid_1Sigmoidwhile/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_1Ќ
while/gru_cell_4/mul_1Mulwhile/gru_cell_4/Sigmoid_1:y:0!while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_1Ј
while/gru_cell_4/add_2AddV2while/gru_cell_4/split:output:2while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_2
while/gru_cell_4/Sigmoid_2Sigmoidwhile/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/Sigmoid_2
while/gru_cell_4/mul_2Mulwhile/gru_cell_4/Sigmoid:y:0while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_2u
while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2
while/gru_cell_4/sub/xЄ
while/gru_cell_4/subSubwhile/gru_cell_4/sub/x:output:0while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/subЃ
while/gru_cell_4/mul_3Mulwhile/gru_cell_4/sub:z:0while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/mul_3Ѓ
while/gru_cell_4/add_3AddV2while/gru_cell_4/mul_2:z:0while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/gru_cell_4/add_3о
*while/TensorArrayV2Write/TensorListSetItemTensorListSetItemwhile_placeholder_1while_placeholderwhile/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02,
*while/TensorArrayV2Write/TensorListSetItem\
while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add/yi
	while/addAddV2while_placeholderwhile/add/y:output:0*
T0*
_output_shapes
: 2
	while/add`
while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2
while/add_1/yv
while/add_1AddV2while_while_loop_counterwhile/add_1/y:output:0*
T0*
_output_shapes
: 2
while/add_1^
while/IdentityIdentitywhile/add_1:z:0*
T0*
_output_shapes
: 2
while/Identityq
while/Identity_1Identitywhile_while_maximum_iterations*
T0*
_output_shapes
: 2
while/Identity_1`
while/Identity_2Identitywhile/add:z:0*
T0*
_output_shapes
: 2
while/Identity_2
while/Identity_3Identity:while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2
while/Identity_3~
while/Identity_4Identitywhile/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2
while/Identity_4"h
1while_gru_cell_4_matmul_1_readvariableop_resource3while_gru_cell_4_matmul_1_readvariableop_resource_0"d
/while_gru_cell_4_matmul_readvariableop_resource1while_gru_cell_4_matmul_readvariableop_resource_0"V
(while_gru_cell_4_readvariableop_resource*while_gru_cell_4_readvariableop_resource_0")
while_identitywhile/Identity:output:0"-
while_identity_1while/Identity_1:output:0"-
while_identity_2while/Identity_2:output:0"-
while_identity_3while/Identity_3:output:0"-
while_identity_4while/Identity_4:output:0"0
while_strided_slice_1while_strided_slice_1_0"Ј
Qwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensorSwhile_tensorarrayv2read_tensorlistgetitem_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
фa
Џ

#sequential_2_gru_4_while_body_67764B
>sequential_2_gru_4_while_sequential_2_gru_4_while_loop_counterH
Dsequential_2_gru_4_while_sequential_2_gru_4_while_maximum_iterations(
$sequential_2_gru_4_while_placeholder*
&sequential_2_gru_4_while_placeholder_1*
&sequential_2_gru_4_while_placeholder_2A
=sequential_2_gru_4_while_sequential_2_gru_4_strided_slice_1_0}
ysequential_2_gru_4_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_4_tensorarrayunstack_tensorlistfromtensor_0A
=sequential_2_gru_4_while_gru_cell_4_readvariableop_resource_0H
Dsequential_2_gru_4_while_gru_cell_4_matmul_readvariableop_resource_0J
Fsequential_2_gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0%
!sequential_2_gru_4_while_identity'
#sequential_2_gru_4_while_identity_1'
#sequential_2_gru_4_while_identity_2'
#sequential_2_gru_4_while_identity_3'
#sequential_2_gru_4_while_identity_4?
;sequential_2_gru_4_while_sequential_2_gru_4_strided_slice_1{
wsequential_2_gru_4_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_4_tensorarrayunstack_tensorlistfromtensor?
;sequential_2_gru_4_while_gru_cell_4_readvariableop_resourceF
Bsequential_2_gru_4_while_gru_cell_4_matmul_readvariableop_resourceH
Dsequential_2_gru_4_while_gru_cell_4_matmul_1_readvariableop_resourceщ
Jsequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2L
Jsequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shapeХ
<sequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItemTensorListGetItemysequential_2_gru_4_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_4_tensorarrayunstack_tensorlistfromtensor_0$sequential_2_gru_4_while_placeholderSsequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItem/element_shape:output:0*'
_output_shapes
:џџџџџџџџџ*
element_dtype02>
<sequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItemн
3sequential_2/gru_4/while/gru_cell_4/ones_like/ShapeShapeCsequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItem:item:0*
T0*
_output_shapes
:25
3sequential_2/gru_4/while/gru_cell_4/ones_like/ShapeЏ
3sequential_2/gru_4/while/gru_cell_4/ones_like/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *  ?25
3sequential_2/gru_4/while/gru_cell_4/ones_like/Const
-sequential_2/gru_4/while/gru_cell_4/ones_likeFill<sequential_2/gru_4/while/gru_cell_4/ones_like/Shape:output:0<sequential_2/gru_4/while/gru_cell_4/ones_like/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_4/while/gru_cell_4/ones_likeц
2sequential_2/gru_4/while/gru_cell_4/ReadVariableOpReadVariableOp=sequential_2_gru_4_while_gru_cell_4_readvariableop_resource_0*
_output_shapes

:0*
dtype024
2sequential_2/gru_4/while/gru_cell_4/ReadVariableOpж
+sequential_2/gru_4/while/gru_cell_4/unstackUnpack:sequential_2/gru_4/while/gru_cell_4/ReadVariableOp:value:0*
T0* 
_output_shapes
:0:0*	
num2-
+sequential_2/gru_4/while/gru_cell_4/unstack
'sequential_2/gru_4/while/gru_cell_4/mulMulCsequential_2/gru_4/while/TensorArrayV2Read/TensorListGetItem:item:06sequential_2/gru_4/while/gru_cell_4/ones_like:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/while/gru_cell_4/mulћ
9sequential_2/gru_4/while/gru_cell_4/MatMul/ReadVariableOpReadVariableOpDsequential_2_gru_4_while_gru_cell_4_matmul_readvariableop_resource_0*
_output_shapes

:0*
dtype02;
9sequential_2/gru_4/while/gru_cell_4/MatMul/ReadVariableOp
*sequential_2/gru_4/while/gru_cell_4/MatMulMatMul+sequential_2/gru_4/while/gru_cell_4/mul:z:0Asequential_2/gru_4/while/gru_cell_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02,
*sequential_2/gru_4/while/gru_cell_4/MatMul
+sequential_2/gru_4/while/gru_cell_4/BiasAddBiasAdd4sequential_2/gru_4/while/gru_cell_4/MatMul:product:04sequential_2/gru_4/while/gru_cell_4/unstack:output:0*
T0*'
_output_shapes
:џџџџџџџџџ02-
+sequential_2/gru_4/while/gru_cell_4/BiasAdd
)sequential_2/gru_4/while/gru_cell_4/ConstConst*
_output_shapes
: *
dtype0*
value	B :2+
)sequential_2/gru_4/while/gru_cell_4/ConstЕ
3sequential_2/gru_4/while/gru_cell_4/split/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ25
3sequential_2/gru_4/while/gru_cell_4/split/split_dimМ
)sequential_2/gru_4/while/gru_cell_4/splitSplit<sequential_2/gru_4/while/gru_cell_4/split/split_dim:output:04sequential_2/gru_4/while/gru_cell_4/BiasAdd:output:0*
T0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2+
)sequential_2/gru_4/while/gru_cell_4/split
;sequential_2/gru_4/while/gru_cell_4/MatMul_1/ReadVariableOpReadVariableOpFsequential_2_gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0*
_output_shapes

:0*
dtype02=
;sequential_2/gru_4/while/gru_cell_4/MatMul_1/ReadVariableOp
,sequential_2/gru_4/while/gru_cell_4/MatMul_1MatMul&sequential_2_gru_4_while_placeholder_2Csequential_2/gru_4/while/gru_cell_4/MatMul_1/ReadVariableOp:value:0*
T0*'
_output_shapes
:џџџџџџџџџ02.
,sequential_2/gru_4/while/gru_cell_4/MatMul_1
-sequential_2/gru_4/while/gru_cell_4/BiasAdd_1BiasAdd6sequential_2/gru_4/while/gru_cell_4/MatMul_1:product:04sequential_2/gru_4/while/gru_cell_4/unstack:output:1*
T0*'
_output_shapes
:џџџџџџџџџ02/
-sequential_2/gru_4/while/gru_cell_4/BiasAdd_1Џ
+sequential_2/gru_4/while/gru_cell_4/Const_1Const*
_output_shapes
:*
dtype0*!
valueB"      џџџџ2-
+sequential_2/gru_4/while/gru_cell_4/Const_1Й
5sequential_2/gru_4/while/gru_cell_4/split_1/split_dimConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ27
5sequential_2/gru_4/while/gru_cell_4/split_1/split_dim
+sequential_2/gru_4/while/gru_cell_4/split_1SplitV6sequential_2/gru_4/while/gru_cell_4/BiasAdd_1:output:04sequential_2/gru_4/while/gru_cell_4/Const_1:output:0>sequential_2/gru_4/while/gru_cell_4/split_1/split_dim:output:0*
T0*

Tlen0*M
_output_shapes;
9:џџџџџџџџџ:џџџџџџџџџ:џџџџџџџџџ*
	num_split2-
+sequential_2/gru_4/while/gru_cell_4/split_1ї
'sequential_2/gru_4/while/gru_cell_4/addAddV22sequential_2/gru_4/while/gru_cell_4/split:output:04sequential_2/gru_4/while/gru_cell_4/split_1:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/while/gru_cell_4/addФ
+sequential_2/gru_4/while/gru_cell_4/SigmoidSigmoid+sequential_2/gru_4/while/gru_cell_4/add:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2-
+sequential_2/gru_4/while/gru_cell_4/Sigmoidћ
)sequential_2/gru_4/while/gru_cell_4/add_1AddV22sequential_2/gru_4/while/gru_cell_4/split:output:14sequential_2/gru_4/while/gru_cell_4/split_1:output:1*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/add_1Ъ
-sequential_2/gru_4/while/gru_cell_4/Sigmoid_1Sigmoid-sequential_2/gru_4/while/gru_cell_4/add_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_4/while/gru_cell_4/Sigmoid_1ј
)sequential_2/gru_4/while/gru_cell_4/mul_1Mul1sequential_2/gru_4/while/gru_cell_4/Sigmoid_1:y:04sequential_2/gru_4/while/gru_cell_4/split_1:output:2*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/mul_1є
)sequential_2/gru_4/while/gru_cell_4/add_2AddV22sequential_2/gru_4/while/gru_cell_4/split:output:2-sequential_2/gru_4/while/gru_cell_4/mul_1:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/add_2Ъ
-sequential_2/gru_4/while/gru_cell_4/Sigmoid_2Sigmoid-sequential_2/gru_4/while/gru_cell_4/add_2:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2/
-sequential_2/gru_4/while/gru_cell_4/Sigmoid_2ш
)sequential_2/gru_4/while/gru_cell_4/mul_2Mul/sequential_2/gru_4/while/gru_cell_4/Sigmoid:y:0&sequential_2_gru_4_while_placeholder_2*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/mul_2
)sequential_2/gru_4/while/gru_cell_4/sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  ?2+
)sequential_2/gru_4/while/gru_cell_4/sub/x№
'sequential_2/gru_4/while/gru_cell_4/subSub2sequential_2/gru_4/while/gru_cell_4/sub/x:output:0/sequential_2/gru_4/while/gru_cell_4/Sigmoid:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2)
'sequential_2/gru_4/while/gru_cell_4/subя
)sequential_2/gru_4/while/gru_cell_4/mul_3Mul+sequential_2/gru_4/while/gru_cell_4/sub:z:01sequential_2/gru_4/while/gru_cell_4/Sigmoid_2:y:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/mul_3я
)sequential_2/gru_4/while/gru_cell_4/add_3AddV2-sequential_2/gru_4/while/gru_cell_4/mul_2:z:0-sequential_2/gru_4/while/gru_cell_4/mul_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2+
)sequential_2/gru_4/while/gru_cell_4/add_3Н
=sequential_2/gru_4/while/TensorArrayV2Write/TensorListSetItemTensorListSetItem&sequential_2_gru_4_while_placeholder_1$sequential_2_gru_4_while_placeholder-sequential_2/gru_4/while/gru_cell_4/add_3:z:0*
_output_shapes
: *
element_dtype02?
=sequential_2/gru_4/while/TensorArrayV2Write/TensorListSetItem
sequential_2/gru_4/while/add/yConst*
_output_shapes
: *
dtype0*
value	B :2 
sequential_2/gru_4/while/add/yЕ
sequential_2/gru_4/while/addAddV2$sequential_2_gru_4_while_placeholder'sequential_2/gru_4/while/add/y:output:0*
T0*
_output_shapes
: 2
sequential_2/gru_4/while/add
 sequential_2/gru_4/while/add_1/yConst*
_output_shapes
: *
dtype0*
value	B :2"
 sequential_2/gru_4/while/add_1/yе
sequential_2/gru_4/while/add_1AddV2>sequential_2_gru_4_while_sequential_2_gru_4_while_loop_counter)sequential_2/gru_4/while/add_1/y:output:0*
T0*
_output_shapes
: 2 
sequential_2/gru_4/while/add_1
!sequential_2/gru_4/while/IdentityIdentity"sequential_2/gru_4/while/add_1:z:0*
T0*
_output_shapes
: 2#
!sequential_2/gru_4/while/IdentityН
#sequential_2/gru_4/while/Identity_1IdentityDsequential_2_gru_4_while_sequential_2_gru_4_while_maximum_iterations*
T0*
_output_shapes
: 2%
#sequential_2/gru_4/while/Identity_1
#sequential_2/gru_4/while/Identity_2Identity sequential_2/gru_4/while/add:z:0*
T0*
_output_shapes
: 2%
#sequential_2/gru_4/while/Identity_2Ц
#sequential_2/gru_4/while/Identity_3IdentityMsequential_2/gru_4/while/TensorArrayV2Write/TensorListSetItem:output_handle:0*
T0*
_output_shapes
: 2%
#sequential_2/gru_4/while/Identity_3З
#sequential_2/gru_4/while/Identity_4Identity-sequential_2/gru_4/while/gru_cell_4/add_3:z:0*
T0*'
_output_shapes
:џџџџџџџџџ2%
#sequential_2/gru_4/while/Identity_4"
Dsequential_2_gru_4_while_gru_cell_4_matmul_1_readvariableop_resourceFsequential_2_gru_4_while_gru_cell_4_matmul_1_readvariableop_resource_0"
Bsequential_2_gru_4_while_gru_cell_4_matmul_readvariableop_resourceDsequential_2_gru_4_while_gru_cell_4_matmul_readvariableop_resource_0"|
;sequential_2_gru_4_while_gru_cell_4_readvariableop_resource=sequential_2_gru_4_while_gru_cell_4_readvariableop_resource_0"O
!sequential_2_gru_4_while_identity*sequential_2/gru_4/while/Identity:output:0"S
#sequential_2_gru_4_while_identity_1,sequential_2/gru_4/while/Identity_1:output:0"S
#sequential_2_gru_4_while_identity_2,sequential_2/gru_4/while/Identity_2:output:0"S
#sequential_2_gru_4_while_identity_3,sequential_2/gru_4/while/Identity_3:output:0"S
#sequential_2_gru_4_while_identity_4,sequential_2/gru_4/while/Identity_4:output:0"|
;sequential_2_gru_4_while_sequential_2_gru_4_strided_slice_1=sequential_2_gru_4_while_sequential_2_gru_4_strided_slice_1_0"є
wsequential_2_gru_4_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_4_tensorarrayunstack_tensorlistfromtensorysequential_2_gru_4_while_tensorarrayv2read_tensorlistgetitem_sequential_2_gru_4_tensorarrayunstack_tensorlistfromtensor_0*>
_input_shapes-
+: : : : :џџџџџџџџџ: : :::: 

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :

_output_shapes
: :-)
'
_output_shapes
:џџџџџџџџџ:

_output_shapes
: :

_output_shapes
: 
Ь<
Ю
@__inference_gru_4_layer_call_and_return_conditional_losses_68495

inputs
gru_cell_4_68419
gru_cell_4_68421
gru_cell_4_68423
identityЂ"gru_cell_4/StatefulPartitionedCallЂwhileD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2т
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice\
zeros/mul/yConst*
_output_shapes
: *
dtype0*
value	B :2
zeros/mul/yl
	zeros/mulMulstrided_slice:output:0zeros/mul/y:output:0*
T0*
_output_shapes
: 2
	zeros/mul_
zeros/Less/yConst*
_output_shapes
: *
dtype0*
value
B :ш2
zeros/Less/yg

zeros/LessLesszeros/mul:z:0zeros/Less/y:output:0*
T0*
_output_shapes
: 2

zeros/Lessb
zeros/packed/1Const*
_output_shapes
: *
dtype0*
value	B :2
zeros/packed/1
zeros/packedPackstrided_slice:output:0zeros/packed/1:output:0*
N*
T0*
_output_shapes
:2
zeros/packed_
zeros/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *    2
zeros/Constu
zerosFillzeros/packed:output:0zeros/Const:output:0*
T0*'
_output_shapes
:џџџџџџџџџ2
zerosu
transpose/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose/perm
	transpose	Transposeinputstranspose/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
	transposeO
Shape_1Shapetranspose:y:0*
T0*
_output_shapes
:2	
Shape_1x
strided_slice_1/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_1/stack|
strided_slice_1/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_1|
strided_slice_1/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_1/stack_2ю
strided_slice_1StridedSliceShape_1:output:0strided_slice_1/stack:output:0 strided_slice_1/stack_1:output:0 strided_slice_1/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_slice_1
TensorArrayV2/element_shapeConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
TensorArrayV2/element_shapeВ
TensorArrayV2TensorListReserve$TensorArrayV2/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2П
5TensorArrayUnstack/TensorListFromTensor/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   27
5TensorArrayUnstack/TensorListFromTensor/element_shapeј
'TensorArrayUnstack/TensorListFromTensorTensorListFromTensortranspose:y:0>TensorArrayUnstack/TensorListFromTensor/element_shape:output:0*
_output_shapes
: *
element_dtype0*

shape_type02)
'TensorArrayUnstack/TensorListFromTensorx
strided_slice_2/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_2/stack|
strided_slice_2/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_1|
strided_slice_2/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_2/stack_2ќ
strided_slice_2StridedSlicetranspose:y:0strided_slice_2/stack:output:0 strided_slice_2/stack_1:output:0 strided_slice_2/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_2ц
"gru_cell_4/StatefulPartitionedCallStatefulPartitionedCallstrided_slice_2:output:0zeros:output:0gru_cell_4_68419gru_cell_4_68421gru_cell_4_68423*
Tin	
2*
Tout
2*
_collective_manager_ids
 *:
_output_shapes(
&:џџџџџџџџџ:џџџџџџџџџ*%
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8 *N
fIRG
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_681282$
"gru_cell_4/StatefulPartitionedCall
TensorArrayV2_1/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   2
TensorArrayV2_1/element_shapeИ
TensorArrayV2_1TensorListReserve&TensorArrayV2_1/element_shape:output:0strided_slice_1:output:0*
_output_shapes
: *
element_dtype0*

shape_type02
TensorArrayV2_1N
timeConst*
_output_shapes
: *
dtype0*
value	B : 2
time
while/maximum_iterationsConst*
_output_shapes
: *
dtype0*
valueB :
џџџџџџџџџ2
while/maximum_iterationsj
while/loop_counterConst*
_output_shapes
: *
dtype0*
value	B : 2
while/loop_counterп
whileWhilewhile/loop_counter:output:0!while/maximum_iterations:output:0time:output:0TensorArrayV2_1:handle:0zeros:output:0strided_slice_1:output:07TensorArrayUnstack/TensorListFromTensor:output_handle:0gru_cell_4_68419gru_cell_4_68421gru_cell_4_68423*
T
2
*
_lower_using_switch_merge(*
_num_original_outputs
*9
_output_shapes'
%: : : : :џџџџџџџџџ: : : : : *%
_read_only_resource_inputs
	*
bodyR
while_body_68431*
condR
while_cond_68430*8
output_shapes'
%: : : : :џџџџџџџџџ: : : : : *
parallel_iterations 2
whileЕ
0TensorArrayV2Stack/TensorListStack/element_shapeConst*
_output_shapes
:*
dtype0*
valueB"џџџџ   22
0TensorArrayV2Stack/TensorListStack/element_shapeё
"TensorArrayV2Stack/TensorListStackTensorListStackwhile:output:39TensorArrayV2Stack/TensorListStack/element_shape:output:0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ*
element_dtype02$
"TensorArrayV2Stack/TensorListStack
strided_slice_3/stackConst*
_output_shapes
:*
dtype0*
valueB:
џџџџџџџџџ2
strided_slice_3/stack|
strided_slice_3/stack_1Const*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice_3/stack_1|
strided_slice_3/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice_3/stack_2
strided_slice_3StridedSlice+TensorArrayV2Stack/TensorListStack:tensor:0strided_slice_3/stack:output:0 strided_slice_3/stack_1:output:0 strided_slice_3/stack_2:output:0*
Index0*
T0*'
_output_shapes
:џџџџџџџџџ*
shrink_axis_mask2
strided_slice_3y
transpose_1/permConst*
_output_shapes
:*
dtype0*!
valueB"          2
transpose_1/permЎ
transpose_1	Transpose+TensorArrayV2Stack/TensorListStack:tensor:0transpose_1/perm:output:0*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2
transpose_1f
runtimeConst"/device:CPU:0*
_output_shapes
: *
dtype0*
valueB
 *    2	
runtime
IdentityIdentitytranspose_1:y:0#^gru_cell_4/StatefulPartitionedCall^while*
T0*4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ2

Identity"
identityIdentity:output:0*?
_input_shapes.
,:џџџџџџџџџџџџџџџџџџ:::2H
"gru_cell_4/StatefulPartitionedCall"gru_cell_4/StatefulPartitionedCall2
whilewhile:\ X
4
_output_shapes"
 :џџџџџџџџџџџџџџџџџџ
 
_user_specified_nameinputs"ИL
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*П
serving_defaultЋ
P
embedding_2_input;
#serving_default_embedding_2_input:0џџџџџџџџџш;
dense_20
StatefulPartitionedCall:0џџџџџџџџџtensorflow/serving/predict:Ё
Є5
layer_with_weights-0
layer-0
layer_with_weights-1
layer-1
layer_with_weights-2
layer-2
layer_with_weights-3
layer-3
	optimizer
trainable_variables
regularization_losses
	variables
		keras_api


signatures
w_default_save_signature
*x&call_and_return_all_conditional_losses
y__call__"Н2
_tf_keras_sequential2{"class_name": "Sequential", "name": "sequential_2", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "must_restore_from_config": false, "config": {"name": "sequential_2", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "embedding_2_input"}}, {"class_name": "Embedding", "config": {"name": "embedding_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "dtype": "float32", "input_dim": 257, "output_dim": 16, "embeddings_initializer": {"class_name": "RandomUniform", "config": {"minval": -0.05, "maxval": 0.05, "seed": null}}, "embeddings_regularizer": null, "activity_regularizer": null, "embeddings_constraint": null, "mask_zero": false, "input_length": 1000}}, {"class_name": "GRU", "config": {"name": "gru_4", "trainable": true, "dtype": "float32", "return_sequences": true, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}, {"class_name": "GRU", "config": {"name": "gru_5", "trainable": true, "dtype": "float32", "return_sequences": false, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1000]}, "is_graph_network": true, "keras_version": "2.4.0", "backend": "tensorflow", "model_config": {"class_name": "Sequential", "config": {"name": "sequential_2", "layers": [{"class_name": "InputLayer", "config": {"batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "dtype": "float32", "sparse": false, "ragged": false, "name": "embedding_2_input"}}, {"class_name": "Embedding", "config": {"name": "embedding_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "dtype": "float32", "input_dim": 257, "output_dim": 16, "embeddings_initializer": {"class_name": "RandomUniform", "config": {"minval": -0.05, "maxval": 0.05, "seed": null}}, "embeddings_regularizer": null, "activity_regularizer": null, "embeddings_constraint": null, "mask_zero": false, "input_length": 1000}}, {"class_name": "GRU", "config": {"name": "gru_4", "trainable": true, "dtype": "float32", "return_sequences": true, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}, {"class_name": "GRU", "config": {"name": "gru_5", "trainable": true, "dtype": "float32", "return_sequences": false, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}, {"class_name": "Dense", "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}}]}}, "training_config": {"loss": {"class_name": "BinaryCrossentropy", "config": {"reduction": "auto", "name": "binary_crossentropy", "from_logits": false, "label_smoothing": 0}}, "metrics": ["accuracy"], "weighted_metrics": null, "loss_weights": null, "optimizer_config": {"class_name": "Adam", "config": {"name": "Adam", "learning_rate": 0.009999999776482582, "decay": 0.0, "beta_1": 0.8999999761581421, "beta_2": 0.9990000128746033, "epsilon": 1e-07, "amsgrad": false}}}}
Џ

embeddings
trainable_variables
regularization_losses
	variables
	keras_api
*z&call_and_return_all_conditional_losses
{__call__"
_tf_keras_layerі{"class_name": "Embedding", "name": "embedding_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "stateful": false, "must_restore_from_config": false, "config": {"name": "embedding_2", "trainable": true, "batch_input_shape": {"class_name": "__tuple__", "items": [null, 1000]}, "dtype": "float32", "input_dim": 257, "output_dim": 16, "embeddings_initializer": {"class_name": "RandomUniform", "config": {"minval": -0.05, "maxval": 0.05, "seed": null}}, "embeddings_regularizer": null, "activity_regularizer": null, "embeddings_constraint": null, "mask_zero": false, "input_length": 1000}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 1000]}}
Й
cell

state_spec
trainable_variables
regularization_losses
	variables
	keras_api
*|&call_and_return_all_conditional_losses
}__call__"

_tf_keras_rnn_layerђ	{"class_name": "GRU", "name": "gru_4", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "gru_4", "trainable": true, "dtype": "float32", "return_sequences": true, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, null, 16]}, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 1000, 16]}}
К
cell

state_spec
trainable_variables
regularization_losses
	variables
	keras_api
*~&call_and_return_all_conditional_losses
__call__"

_tf_keras_rnn_layerѓ	{"class_name": "GRU", "name": "gru_5", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "gru_5", "trainable": true, "dtype": "float32", "return_sequences": false, "return_state": false, "go_backwards": false, "stateful": false, "unroll": false, "time_major": false, "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}, "input_spec": [{"class_name": "InputSpec", "config": {"dtype": null, "shape": {"class_name": "__tuple__", "items": [null, null, 16]}, "ndim": 3, "max_ndim": null, "min_ndim": null, "axes": {}}}], "build_input_shape": {"class_name": "TensorShape", "items": [null, 1000, 16]}}
є

kernel
bias
trainable_variables
regularization_losses
 	variables
!	keras_api
+&call_and_return_all_conditional_losses
__call__"Э
_tf_keras_layerГ{"class_name": "Dense", "name": "dense_2", "trainable": true, "expects_training_arg": false, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "dense_2", "trainable": true, "dtype": "float32", "units": 1, "activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "bias_regularizer": null, "activity_regularizer": null, "kernel_constraint": null, "bias_constraint": null}, "input_spec": {"class_name": "InputSpec", "config": {"dtype": null, "shape": null, "ndim": null, "max_ndim": null, "min_ndim": 2, "axes": {"-1": 16}}}, "build_input_shape": {"class_name": "TensorShape", "items": [null, 16]}}
ѕ
"iter

#beta_1

$beta_2
	%decay
&learning_ratememfmg'mh(mi)mj*mk+ml,mmvnvovp'vq(vr)vs*vt+vu,vv"
	optimizer
_
0
'1
(2
)3
*4
+5
,6
7
8"
trackable_list_wrapper
 "
trackable_list_wrapper
_
0
'1
(2
)3
*4
+5
,6
7
8"
trackable_list_wrapper
Ъ
trainable_variables
-metrics

.layers
regularization_losses
	variables
/layer_regularization_losses
0non_trainable_variables
1layer_metrics
y__call__
w_default_save_signature
*x&call_and_return_all_conditional_losses
&x"call_and_return_conditional_losses"
_generic_user_object
-
serving_default"
signature_map
):'	2embedding_2/embeddings
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
­
trainable_variables
2metrics

3layers
regularization_losses
	variables
4layer_regularization_losses
5non_trainable_variables
6layer_metrics
{__call__
*z&call_and_return_all_conditional_losses
&z"call_and_return_conditional_losses"
_generic_user_object
Ѕ

'kernel
(recurrent_kernel
)bias
7trainable_variables
8regularization_losses
9	variables
:	keras_api
+&call_and_return_all_conditional_losses
__call__"ш
_tf_keras_layerЮ{"class_name": "GRUCell", "name": "gru_cell_4", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "gru_cell_4", "trainable": true, "dtype": "float32", "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}
 "
trackable_list_wrapper
5
'0
(1
)2"
trackable_list_wrapper
 "
trackable_list_wrapper
5
'0
(1
)2"
trackable_list_wrapper
Й
trainable_variables
;metrics

<states

=layers
regularization_losses
	variables
>layer_regularization_losses
?non_trainable_variables
@layer_metrics
}__call__
*|&call_and_return_all_conditional_losses
&|"call_and_return_conditional_losses"
_generic_user_object
Ѕ

*kernel
+recurrent_kernel
,bias
Atrainable_variables
Bregularization_losses
C	variables
D	keras_api
+&call_and_return_all_conditional_losses
__call__"ш
_tf_keras_layerЮ{"class_name": "GRUCell", "name": "gru_cell_5", "trainable": true, "expects_training_arg": true, "dtype": "float32", "batch_input_shape": null, "stateful": false, "must_restore_from_config": false, "config": {"name": "gru_cell_5", "trainable": true, "dtype": "float32", "units": 16, "activation": "sigmoid", "recurrent_activation": "sigmoid", "use_bias": true, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": null}}, "recurrent_initializer": {"class_name": "Orthogonal", "config": {"gain": 1.0, "seed": null}}, "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": null, "recurrent_regularizer": null, "bias_regularizer": null, "kernel_constraint": null, "recurrent_constraint": null, "bias_constraint": null, "dropout": 0.5, "recurrent_dropout": 0.0, "implementation": 2, "reset_after": true}}
 "
trackable_list_wrapper
5
*0
+1
,2"
trackable_list_wrapper
 "
trackable_list_wrapper
5
*0
+1
,2"
trackable_list_wrapper
Й
trainable_variables
Emetrics

Fstates

Glayers
regularization_losses
	variables
Hlayer_regularization_losses
Inon_trainable_variables
Jlayer_metrics
__call__
*~&call_and_return_all_conditional_losses
&~"call_and_return_conditional_losses"
_generic_user_object
 :2dense_2/kernel
:2dense_2/bias
.
0
1"
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
А
trainable_variables
Kmetrics

Llayers
regularization_losses
 	variables
Mlayer_regularization_losses
Nnon_trainable_variables
Olayer_metrics
__call__
+&call_and_return_all_conditional_losses
'"call_and_return_conditional_losses"
_generic_user_object
:	 (2	Adam/iter
: (2Adam/beta_1
: (2Adam/beta_2
: (2
Adam/decay
: (2Adam/learning_rate
):'02gru_4/gru_cell_4/kernel
3:102!gru_4/gru_cell_4/recurrent_kernel
':%02gru_4/gru_cell_4/bias
):'02gru_5/gru_cell_5/kernel
3:102!gru_5/gru_cell_5/recurrent_kernel
':%02gru_5/gru_cell_5/bias
.
P0
Q1"
trackable_list_wrapper
<
0
1
2
3"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
5
'0
(1
)2"
trackable_list_wrapper
 "
trackable_list_wrapper
5
'0
(1
)2"
trackable_list_wrapper
А
7trainable_variables
Rmetrics

Slayers
8regularization_losses
9	variables
Tlayer_regularization_losses
Unon_trainable_variables
Vlayer_metrics
__call__
+&call_and_return_all_conditional_losses
'"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
5
*0
+1
,2"
trackable_list_wrapper
 "
trackable_list_wrapper
5
*0
+1
,2"
trackable_list_wrapper
А
Atrainable_variables
Wmetrics

Xlayers
Bregularization_losses
C	variables
Ylayer_regularization_losses
Znon_trainable_variables
[layer_metrics
__call__
+&call_and_return_all_conditional_losses
'"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
Л
	\total
	]count
^	variables
_	keras_api"
_tf_keras_metricj{"class_name": "Mean", "name": "loss", "dtype": "float32", "config": {"name": "loss", "dtype": "float32"}}
њ
	`total
	acount
b
_fn_kwargs
c	variables
d	keras_api"Г
_tf_keras_metric{"class_name": "MeanMetricWrapper", "name": "accuracy", "dtype": "float32", "config": {"name": "accuracy", "dtype": "float32", "fn": "binary_accuracy"}}
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
:  (2total
:  (2count
.
\0
]1"
trackable_list_wrapper
-
^	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper
.
`0
a1"
trackable_list_wrapper
-
c	variables"
_generic_user_object
.:,	2Adam/embedding_2/embeddings/m
%:#2Adam/dense_2/kernel/m
:2Adam/dense_2/bias/m
.:,02Adam/gru_4/gru_cell_4/kernel/m
8:602(Adam/gru_4/gru_cell_4/recurrent_kernel/m
,:*02Adam/gru_4/gru_cell_4/bias/m
.:,02Adam/gru_5/gru_cell_5/kernel/m
8:602(Adam/gru_5/gru_cell_5/recurrent_kernel/m
,:*02Adam/gru_5/gru_cell_5/bias/m
.:,	2Adam/embedding_2/embeddings/v
%:#2Adam/dense_2/kernel/v
:2Adam/dense_2/bias/v
.:,02Adam/gru_4/gru_cell_4/kernel/v
8:602(Adam/gru_4/gru_cell_4/recurrent_kernel/v
,:*02Adam/gru_4/gru_cell_4/bias/v
.:,02Adam/gru_5/gru_cell_5/kernel/v
8:602(Adam/gru_5/gru_cell_5/recurrent_kernel/v
,:*02Adam/gru_5/gru_cell_5/bias/v
щ2ц
 __inference__wrapped_model_68028С
В
FullArgSpec
args 
varargsjargs
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *1Ђ.
,)
embedding_2_inputџџџџџџџџџш
ъ2ч
G__inference_sequential_2_layer_call_and_return_conditional_losses_70092
G__inference_sequential_2_layer_call_and_return_conditional_losses_70118
G__inference_sequential_2_layer_call_and_return_conditional_losses_70689
G__inference_sequential_2_layer_call_and_return_conditional_losses_71032Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ў2ћ
,__inference_sequential_2_layer_call_fn_70217
,__inference_sequential_2_layer_call_fn_71078
,__inference_sequential_2_layer_call_fn_71055
,__inference_sequential_2_layer_call_fn_70168Р
ЗВГ
FullArgSpec1
args)&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
№2э
F__inference_embedding_2_layer_call_and_return_conditional_losses_71088Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
е2в
+__inference_embedding_2_layer_call_fn_71095Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
у2р
@__inference_gru_4_layer_call_and_return_conditional_losses_71714
@__inference_gru_4_layer_call_and_return_conditional_losses_71477
@__inference_gru_4_layer_call_and_return_conditional_losses_71310
@__inference_gru_4_layer_call_and_return_conditional_losses_71881е
ЬВШ
FullArgSpecB
args:7
jself
jinputs
jmask

jtraining
jinitial_state
varargs
 
varkw
 
defaults

 
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ї2є
%__inference_gru_4_layer_call_fn_71488
%__inference_gru_4_layer_call_fn_71903
%__inference_gru_4_layer_call_fn_71892
%__inference_gru_4_layer_call_fn_71499е
ЬВШ
FullArgSpecB
args:7
jself
jinputs
jmask

jtraining
jinitial_state
varargs
 
varkw
 
defaults

 
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
у2р
@__inference_gru_5_layer_call_and_return_conditional_losses_72118
@__inference_gru_5_layer_call_and_return_conditional_losses_72689
@__inference_gru_5_layer_call_and_return_conditional_losses_72522
@__inference_gru_5_layer_call_and_return_conditional_losses_72285е
ЬВШ
FullArgSpecB
args:7
jself
jinputs
jmask

jtraining
jinitial_state
varargs
 
varkw
 
defaults

 
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ї2є
%__inference_gru_5_layer_call_fn_72307
%__inference_gru_5_layer_call_fn_72711
%__inference_gru_5_layer_call_fn_72296
%__inference_gru_5_layer_call_fn_72700е
ЬВШ
FullArgSpecB
args:7
jself
jinputs
jmask

jtraining
jinitial_state
varargs
 
varkw
 
defaults

 
p 

 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
ь2щ
B__inference_dense_2_layer_call_and_return_conditional_losses_72722Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
б2Ю
'__inference_dense_2_layer_call_fn_72731Ђ
В
FullArgSpec
args
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs 
kwonlydefaults
 
annotationsЊ *
 
<B:
#__inference_signature_wrapper_70250embedding_2_input
в2Я
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72843
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72799О
ЕВБ
FullArgSpec3
args+(
jself
jinputs
jstates

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
2
*__inference_gru_cell_4_layer_call_fn_72871
*__inference_gru_cell_4_layer_call_fn_72857О
ЕВБ
FullArgSpec3
args+(
jself
jinputs
jstates

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
в2Я
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72939
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72983О
ЕВБ
FullArgSpec3
args+(
jself
jinputs
jstates

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
2
*__inference_gru_cell_5_layer_call_fn_73011
*__inference_gru_cell_5_layer_call_fn_72997О
ЕВБ
FullArgSpec3
args+(
jself
jinputs
jstates

jtraining
varargs
 
varkw
 
defaults
p 

kwonlyargs 
kwonlydefaultsЊ 
annotationsЊ *
 
 __inference__wrapped_model_68028{	)'(,*+;Ђ8
1Ђ.
,)
embedding_2_inputџџџџџџџџџш
Њ "1Њ.
,
dense_2!
dense_2џџџџџџџџџЂ
B__inference_dense_2_layer_call_and_return_conditional_losses_72722\/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "%Ђ"

0џџџџџџџџџ
 z
'__inference_dense_2_layer_call_fn_72731O/Ђ,
%Ђ"
 
inputsџџџџџџџџџ
Њ "џџџџџџџџџЋ
F__inference_embedding_2_layer_call_and_return_conditional_losses_71088a0Ђ-
&Ђ#
!
inputsџџџџџџџџџш
Њ "*Ђ'
 
0џџџџџџџџџш
 
+__inference_embedding_2_layer_call_fn_71095T0Ђ-
&Ђ#
!
inputsџџџџџџџџџш
Њ "џџџџџџџџџшЗ
@__inference_gru_4_layer_call_and_return_conditional_losses_71310s)'(@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p

 
Њ "*Ђ'
 
0џџџџџџџџџш
 З
@__inference_gru_4_layer_call_and_return_conditional_losses_71477s)'(@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p 

 
Њ "*Ђ'
 
0џџџџџџџџџш
 Я
@__inference_gru_4_layer_call_and_return_conditional_losses_71714)'(OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p

 
Њ "2Ђ/
(%
0џџџџџџџџџџџџџџџџџџ
 Я
@__inference_gru_4_layer_call_and_return_conditional_losses_71881)'(OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p 

 
Њ "2Ђ/
(%
0џџџџџџџџџџџџџџџџџџ
 
%__inference_gru_4_layer_call_fn_71488f)'(@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p

 
Њ "џџџџџџџџџш
%__inference_gru_4_layer_call_fn_71499f)'(@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p 

 
Њ "џџџџџџџџџшІ
%__inference_gru_4_layer_call_fn_71892})'(OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p

 
Њ "%"џџџџџџџџџџџџџџџџџџІ
%__inference_gru_4_layer_call_fn_71903})'(OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p 

 
Њ "%"џџџџџџџџџџџџџџџџџџС
@__inference_gru_5_layer_call_and_return_conditional_losses_72118},*+OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p

 
Њ "%Ђ"

0џџџџџџџџџ
 С
@__inference_gru_5_layer_call_and_return_conditional_losses_72285},*+OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 В
@__inference_gru_5_layer_call_and_return_conditional_losses_72522n,*+@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p

 
Њ "%Ђ"

0џџџџџџџџџ
 В
@__inference_gru_5_layer_call_and_return_conditional_losses_72689n,*+@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 
%__inference_gru_5_layer_call_fn_72296p,*+OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p

 
Њ "џџџџџџџџџ
%__inference_gru_5_layer_call_fn_72307p,*+OЂL
EЂB
41
/,
inputs/0џџџџџџџџџџџџџџџџџџ

 
p 

 
Њ "џџџџџџџџџ
%__inference_gru_5_layer_call_fn_72700a,*+@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p

 
Њ "џџџџџџџџџ
%__inference_gru_5_layer_call_fn_72711a,*+@Ђ=
6Ђ3
%"
inputsџџџџџџџџџш

 
p 

 
Њ "џџџџџџџџџ
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72799З)'(\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p
Њ "RЂO
HЂE

0/0џџџџџџџџџ
$!

0/1/0џџџџџџџџџ
 
E__inference_gru_cell_4_layer_call_and_return_conditional_losses_72843З)'(\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p 
Њ "RЂO
HЂE

0/0џџџџџџџџџ
$!

0/1/0џџџџџџџџџ
 и
*__inference_gru_cell_4_layer_call_fn_72857Љ)'(\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p
Њ "DЂA

0џџџџџџџџџ
"

1/0џџџџџџџџџи
*__inference_gru_cell_4_layer_call_fn_72871Љ)'(\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p 
Њ "DЂA

0џџџџџџџџџ
"

1/0џџџџџџџџџ
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72939З,*+\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p
Њ "RЂO
HЂE

0/0џџџџџџџџџ
$!

0/1/0џџџџџџџџџ
 
E__inference_gru_cell_5_layer_call_and_return_conditional_losses_72983З,*+\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p 
Њ "RЂO
HЂE

0/0џџџџџџџџџ
$!

0/1/0џџџџџџџџџ
 и
*__inference_gru_cell_5_layer_call_fn_72997Љ,*+\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p
Њ "DЂA

0џџџџџџџџџ
"

1/0џџџџџџџџџи
*__inference_gru_cell_5_layer_call_fn_73011Љ,*+\ЂY
RЂO
 
inputsџџџџџџџџџ
'Ђ$
"
states/0џџџџџџџџџ
p 
Њ "DЂA

0џџџџџџџџџ
"

1/0џџџџџџџџџТ
G__inference_sequential_2_layer_call_and_return_conditional_losses_70092w	)'(,*+CЂ@
9Ђ6
,)
embedding_2_inputџџџџџџџџџш
p

 
Њ "%Ђ"

0џџџџџџџџџ
 Т
G__inference_sequential_2_layer_call_and_return_conditional_losses_70118w	)'(,*+CЂ@
9Ђ6
,)
embedding_2_inputџџџџџџџџџш
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 З
G__inference_sequential_2_layer_call_and_return_conditional_losses_70689l	)'(,*+8Ђ5
.Ђ+
!
inputsџџџџџџџџџш
p

 
Њ "%Ђ"

0џџџџџџџџџ
 З
G__inference_sequential_2_layer_call_and_return_conditional_losses_71032l	)'(,*+8Ђ5
.Ђ+
!
inputsџџџџџџџџџш
p 

 
Њ "%Ђ"

0џџџџџџџџџ
 
,__inference_sequential_2_layer_call_fn_70168j	)'(,*+CЂ@
9Ђ6
,)
embedding_2_inputџџџџџџџџџш
p

 
Њ "џџџџџџџџџ
,__inference_sequential_2_layer_call_fn_70217j	)'(,*+CЂ@
9Ђ6
,)
embedding_2_inputџџџџџџџџџш
p 

 
Њ "џџџџџџџџџ
,__inference_sequential_2_layer_call_fn_71055_	)'(,*+8Ђ5
.Ђ+
!
inputsџџџџџџџџџш
p

 
Њ "џџџџџџџџџ
,__inference_sequential_2_layer_call_fn_71078_	)'(,*+8Ђ5
.Ђ+
!
inputsџџџџџџџџџш
p 

 
Њ "џџџџџџџџџИ
#__inference_signature_wrapper_70250	)'(,*+PЂM
Ђ 
FЊC
A
embedding_2_input,)
embedding_2_inputџџџџџџџџџш"1Њ.
,
dense_2!
dense_2џџџџџџџџџ