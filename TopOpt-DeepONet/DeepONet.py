import numpy as np
import tensorflow as tf

class DeepONet():
    def __init__(self,layers_Trunk,layers_Branch,actfunc):
        self.dtype = "float32"
        tf.keras.backend.set_floatx(self.dtype)
        self.init_NN(layers_Trunk, layers_Branch,actfunc)
        return

    #################    Load Data    ###############    
    def load_data(self):
        self.Y = tf.convert_to_tensor(np.load('Y.npy'), dtype=self.dtype)
        self.Y = np.transpose(self.Y)
        self.U_data = tf.convert_to_tensor(np.load('U.npy'),dtype=self.dtype)

        self.rho_data = tf.convert_to_tensor(np.load('rho.npy'), dtype=self.dtype)

        self.num_points = self.rho_data.shape[0]
        self.num_data = self.rho_data.shape[1]
        return

    #################    Train Network    ###############    
    @tf.function
    def train(self, max_epoch, e_print):
        loss_list = tf.TensorArray(dtype=self.dtype, size=max_epoch, dynamic_size=False)
        
        for e in tf.range(max_epoch):
            (loss) = self.fit_epoch()
            loss_list = loss_list.write(e,loss)
            if tf.math.mod(e,e_print)==0:
                formatted_tensor = tf.strings.format("e: {}\t  Loss: {}",
                                                     (e, loss))
                tf.print(formatted_tensor)
        return (loss_list.stack())
    
    #################    Network initializaion    ##############    
    def init_NN(self,layers_Trunk, layers_Branch,actfunc):
        # Branch net
        inputs_Branch = tf.keras.layers.Input(shape=(2,))
        Branch = inputs_Branch
        for width in layers_Branch[1:-1]:
            Branch = tf.keras.layers.Dense(width, activation=actfunc,
                kernel_initializer='glorot_normal')(Branch)
        
        Branch = tf.keras.layers.Dense(layers_Branch[-1], activation=None,
            kernel_initializer='glorot_normal')(Branch)
        self.Branch = tf.keras.models.Model(inputs_Branch, Branch)
        self.Branch.compiled_metrics == None    

        # Trunk net
        inputs_Trunk = tf.keras.layers.Input(shape=(2,))
        Trunk = inputs_Trunk
        for width in layers_Trunk[1:-1]:
            Trunk = tf.keras.layers.Dense(width, activation=actfunc,
                kernel_initializer='glorot_normal')(Trunk)
      
        Trunk = tf.keras.layers.Dense(layers_Trunk[-1], activation=None,
            kernel_initializer='glorot_normal')(Trunk)
        self.Trunk = tf.keras.models.Model(inputs_Trunk, Trunk)
        self.Trunk.compiled_metrics == None    
        return

    #################    Optimizer Initialization    ##############        
    def init_optimizer(self, lr, b1):
        self.tf_optimizer = tf.keras.optimizers.Adam(
            learning_rate=lr, beta_1=b1, epsilon = 1e-7)
        return

    #################    Neural network optimization (Fitting)    ######
    def fit_epoch(self):         
        with tf.GradientTape() as tape:
            loss = self.loss_DeepONet()        
            grads = tape.gradient(loss, self.Branch.trainable_variables+
                                        self.Trunk.trainable_variables)
        del tape            
        # Update
        self.tf_optimizer.apply_gradients(
              zip(grads,self.Branch.trainable_variables+
                        self.Trunk.trainable_variables
                 )
        )     
        return loss

    #################    DeepONet Loss Calculation   ######################
    def loss_DeepONet(self):
        loss_DeepONet = 0 
        for i in range(self.num_data):
            rho_data_N = tf.expand_dims(tf.gather(self.rho_data, i, axis=1),axis=0)                
            Y_i = self.Y
            U_i = tf.expand_dims(tf.gather(self.U_data, i, axis=1),axis=0)        
            rho_NN = self.cal_rho_N(U_i,Y_i)
            loss_DeepONet = loss_DeepONet + (tf.reduce_mean(tf.abs(rho_data_N-rho_NN))) / self.num_data
        return loss_DeepONet
                
    #################    NN Field definition    ######################
    def cal_rho_N(self,U_i,Y_i):
        branch = self.Branch(U_i)
        trunk = self.Trunk(Y_i)

        G_1 = tf.math.multiply(branch,trunk)
        rho = tf.reduce_sum(G_1, axis=1)   # Just 1 case
        return rho
