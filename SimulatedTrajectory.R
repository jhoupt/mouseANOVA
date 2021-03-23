minimum_jerk_policy <- function(initial_state, target_state, D, 
				x_tol, y_tol, max_acceleration) { 
  
  qdot <- vector("list")
  Q <- matrix(c(0,1,0,0,0,1,-60/D^3,-36/D^2,-9/D), 
      	3, 3, byrow=TRUE)
  df <- matrix(c(0,0,60/D^3), nrow=3)
  qdot$x <- Q %*% initial_state$x + df * target_state$x
  qdot$y <- Q %*% initial_state$y + df * target_state$y


   
  halt_policy <- (abs(initial_state$x[1] - target_state$x) < x_tol & 
                  abs(initial_state$y[1] - target_state$y) < y_tol & 
                  abs(initial_state$x[2]) < max_acceleration/2 & 
                  abs(initial_state$y[2]) < max_acceleration/2)
  

  if (!is.logical(halt_policy)) { 
    if(any(is.nan(qdot$x)) | any(is.nan(qdot$y))) {
      print(paste("Halt", halt_policy, sep=":"))
    }
    #print(initial_state)
    halt_policy <- FALSE
  }

  tryCatch({
    if (halt_policy) { 
      qdot$x[3] <- -qdot$x[2]
      qdot$y[3] <- -qdot$y[2]
      qdot$x[2] <- -qdot$x[1]
      qdot$y[2] <- -qdot$y[1]
    }
  }, error = function(e) { print(">>>"); print(initial_state)})

  # Recover from error in calculating acceleration
  if(is.nan(qdot$x[3])) {
    qdot$x[3] <- 0
  }
  if(is.nan(qdot$y[3])) {
    qdot$y[3] <- 0
  }
  return(qdot)
}


test_jerk_policy <- function(initial_state, final_state, D, 
			     x_tol, y_tol, max_acceleration) { 
  qdot <- vector("list")
  state <- vector('list', length=D)
  good <- TRUE

  tf <- D
  state[[1]] <- initial_state
  for (tx in 2:tf) { 
    D <- tf - tx
    qdot <- minimum_jerk_policy(state[[tx-1]], final_state, D, x_tol, y_tol,
			        max_acceleration)
    state[[tx]] <- vector("list")
    state[[tx]]$x <- state[[tx-1]]$x + qdot$x
    state[[tx]]$y <- state[[tx-1]]$y + qdot$y
    acc <- sqrt(state[[tx]]$x[3]^2 + state[[tx]]$y[3]^2)
    if (acc > max_acceleration) {
      return(FALSE)
    }
  }
  return(TRUE)
}


rtrajectory <- function(n, initial_acceleration, max_acceleration, 
                        initial_pos, final_pos, sampling_freq, max_t,
                        control_noise, first_deviation) { 
  x_tol = 20
  y_tol = 10 
  
  # Convert to unit steps
  max_acceleration <- max_acceleration * (1/sampling_freq)^2
  initial_acceleration <- initial_acceleration * (1/sampling_freq)^2
  first_deviation <- first_deviation * sampling_freq
  max_t <- max_t * sampling_freq

  sampling_freq <- 1
  delta_t = 1/sampling_freq
  delta_t2 = delta_t^2
  t_vec = seq(0,max_t,by=1/sampling_freq)

  final_state <- vector("list")
  final_state$x <- final_pos[1]
  final_state$y <- final_pos[2]
  
  traj_mat <- array(0, c(n, 2, length(t_vec)))
  vel_mat <- array(0, c(n, 22, length(t_vec)))
  
  for (n in 1:n) { 
    # Initial trajectory is straight upward for a random amount of time
    initial_deviation <- runif(1, first_deviation/2, first_deviation * 1.5)
  
    y_state <- rep("", length(t_vec))
    x_state <- rep("", length(t_vec))
    
    theta = rep(0, length(t_vec))
    initial_steps <- (1+floor(sampling_freq*initial_deviation))

    state <- vector('list', length=length(t_vec))
    
    
    ## Initial indecisive movement ##
    x_state[1] ="Hold"
    y_state[1] = "Hold"
    state[[1]]$x <- matrix(rep(0,3), nrow=3)
    state[[1]]$y <- matrix(rep(0,3), nrow=3)
    tx = 2

    for (tx in 2:initial_steps) { 
      x_state[tx] <- "Hold"
      state[[tx]]$x <- matrix(rep(0,3), nrow=3)
      state[[tx]]$x[3] <- rnorm(1, 0, sd = control_noise)
      state[[tx]]$x[2] <- state[[tx-1]]$x[3] + state[[tx-1]]$x[2]
      state[[tx]]$x[1] <- .5 * state[[tx-1]]$x[3] + state[[tx-1]]$x[2] + 
      	                  state[[tx-1]]$x[1]

      y_state[tx] <- "Accelerate"
      state[[tx]]$y <- matrix(rep(0,3), nrow=3)
      state[[tx]]$y[3] <- initial_acceleration + 
	                  rnorm(1, 0, sd = control_noise)
      state[[tx]]$y[2] <- state[[tx-1]]$y[3] + state[[tx-1]]$y[2]
      state[[tx]]$y[1] <- .5 * state[[tx-1]]$y[3] + state[[tx-1]]$y[2] + 
      	                  state[[tx-1]]$y[1]
    }


    tx <- tx + 1 
    # Hack so that transition from y-accerlation is not counted in jerk
    tmp_x <- state[[tx-1]]$x[3]
    tmp_y <- state[[tx-1]]$y[3]
    state[[tx-1]]$x[3] <- 0 
    state[[tx-1]]$y[3] <- 0

    D <- Inf
    while (tx <= length(t_vec)) {
      if (D > 2) { 
        D <- 3
        check <- FALSE
        while(check == FALSE & D < max_t/delta_t -tx) { 
          check <- test_jerk_policy(state[[tx-1]], final_state, D, 
        		            x_tol, y_tol, max_acceleration)
          D <- D + 1
        }
      }
      if (D > 0) { D <- D-1 }  
      #if (is.nan(state[[tx-1]]$x)) { 
      #  print(state[[tx-3]]$x)
      #  print(state[[tx-2]]$x)
      #  print(state[[tx-3]]$y)
      #  print(state[[tx-2]]$y)
      #}
      control <- minimum_jerk_policy(state[[tx-1]], final_state, D=D,
      			       x_tol, y_tol, max_acceleration)
      state[[tx]] <- vector("list")

      if(any(is.nan(control$x)) | any(is.nan(control$y))) {
        state[[tx]]$x <- matrix(rep(NaN, 3), nrow=3)
        state[[tx]]$y <- matrix(rep(NaN, 3), nrow=3)
      }
      #cat(paste(D), ".")

      state[[tx]]$x <- state[[tx-1]]$x + control$x
      state[[tx]]$y <- state[[tx-1]]$y + control$y

      # Add motor noise
      state[[tx]]$x[3] <- state[[tx]]$x[3] + 
      	            rnorm(1, mean=0, sd=control_noise)
      state[[tx]]$y[3] <- state[[tx]]$y[3] + 
      	            rnorm(1, mean=0, sd=control_noise)
      tx <- tx + 1
    }
    #cat("\n")

    acc_x <- rep(NA, length(t_vec))
    vel_x <- rep(NA, length(t_vec))
    pos_x <- rep(NA, length(t_vec))
    acc_y <- rep(NA, length(t_vec))
    vel_y <- rep(NA, length(t_vec))
    pos_y <- rep(NA, length(t_vec))
    for( tx in 1:length(t_vec)){ 
      pos_x[tx] <- state[[tx]]$x[1]
      vel_x[tx] <- state[[tx]]$x[2]
      acc_x[tx] <- state[[tx]]$x[3]
      pos_y[tx] <- state[[tx]]$y[1]
      vel_y[tx] <- state[[tx]]$y[2]
      acc_y[tx] <- state[[tx]]$y[3]
    }
  
    traj_mat[n,1,] <- pos_x
    traj_mat[n,2,] <- pos_y
  
    vel_mat[n,1,] <- vel_x
    vel_mat[n,2,] <- vel_y
  }
  return(traj_mat)
} 

rtrajectory_decnoise <- function(n, initial_acceleration, max_acceleration, 
                        initial_pos, final_pos, sampling_freq, max_t,
                        control_noise, dec_noise, first_deviation) { 
  x_tol = 20
  y_tol = 10 
  
  # Convert to unit steps
  max_acceleration <- max_acceleration * (1/sampling_freq)^2
  initial_acceleration <- initial_acceleration * (1/sampling_freq)^2
  first_deviation <- first_deviation * sampling_freq
  max_t <- max_t * sampling_freq

  sampling_freq <- 1
  delta_t = 1/sampling_freq
  delta_t2 = delta_t^2
  t_vec = seq(0,max_t,by=1/sampling_freq)

  final_state <- vector("list")
  final_state$x <- final_pos[1]
  final_state$y <- final_pos[2]
  
  traj_mat <- array(0, c(n, 2, length(t_vec)))
  vel_mat <- array(0, c(n, 22, length(t_vec)))
  
  for (n in 1:n) { 
    # Initial trajectory is straight upward for a random amount of time
    initial_deviation <- runif(1, first_deviation/2, first_deviation * 1.5)
  
    y_state <- rep("", length(t_vec))
    x_state <- rep("", length(t_vec))
    
    theta = rep(0, length(t_vec))
    initial_steps <- (1+floor(sampling_freq*initial_deviation))

    state <- vector('list', length=length(t_vec))
    
    
    ## Initial indecisive movement ##
    x_state[1] ="Hold"
    y_state[1] = "Hold"
    state[[1]]$x <- matrix(rep(0,3), nrow=3)
    state[[1]]$y <- matrix(rep(0,3), nrow=3)
    tx = 2

    #for (tx in 2:initial_steps) { 
    #  x_state[tx] <- "Hold"
    #  state[[tx]]$x <- matrix(rep(0,3), nrow=3)
    #  state[[tx]]$x[3] <- rnorm(1, 0, sd = control_noise)
    #  state[[tx]]$x[2] <- state[[tx-1]]$x[3] + state[[tx-1]]$x[2]
    #  state[[tx]]$x[1] <- .5 * state[[tx-1]]$x[3] + state[[tx-1]]$x[2] + 
    #  	                  state[[tx-1]]$x[1]

    #  y_state[tx] <- "Accelerate"
    #  state[[tx]]$y <- matrix(rep(0,3), nrow=3)
    #  state[[tx]]$y[3] <- initial_acceleration + 
    #                      rnorm(1, 0, sd = control_noise)
    #  state[[tx]]$y[2] <- state[[tx-1]]$y[3] + state[[tx-1]]$y[2]
    #  state[[tx]]$y[1] <- .5 * state[[tx-1]]$y[3] + state[[tx-1]]$y[2] + 
    #  	                  state[[tx-1]]$y[1]
    #}


    #tx <- tx + 1 
    ## Hack so that transition from y-acceleration is not counted in jerk
    #tmp_x <- state[[tx-1]]$x[3]
    #tmp_y <- state[[tx-1]]$y[3]
    #state[[tx-1]]$x[3] <- 0 
    #state[[tx-1]]$y[3] <- 0

    D <- Inf
    target_state = list(x=0, y=final_state$y)

    while (tx <= length(t_vec)) {
      if (tx < initial_steps) { 
        if (rpois(1,dec_noise)) {
          if (target_state$x ==0) { 
            target_state$x = final_state$x
            if( runif(1) < .5) {
              target_state$x = -target_state$x
            }
          } else { 
            target_state$x = -target_state$x
          }
        }
      } else if(tx==initial_steps) { 
        target_state$x = final_state$x
      }
      if (D > 2) { 
        D <- 3
        check <- FALSE
        while(check == FALSE & D < max_t/delta_t -tx) { 
          check <- test_jerk_policy(state[[tx-1]], target_state, D, 
        		            x_tol, y_tol, max_acceleration)
          D <- D + 1
        }
      }
      if (D > 0) { D <- D-1 }  
      #if (is.nan(state[[tx-1]]$x)) { 
      #  print(state[[tx-3]]$x)
      #  print(state[[tx-2]]$x)
      #  print(state[[tx-3]]$y)
      #  print(state[[tx-2]]$y)
      #}
      control <- minimum_jerk_policy(state[[tx-1]], target_state, D=D,
      			       x_tol, y_tol, max_acceleration)
      state[[tx]] <- vector("list")

      if(any(is.nan(control$x)) | any(is.nan(control$y))) {
        state[[tx]]$x <- matrix(rep(NaN, 3), nrow=3)
        state[[tx]]$y <- matrix(rep(NaN, 3), nrow=3)
      }
      #cat(paste(D), ".")

      state[[tx]]$x <- state[[tx-1]]$x + control$x
      state[[tx]]$y <- state[[tx-1]]$y + control$y

      # Add motor noise
      state[[tx]]$x[3] <- state[[tx]]$x[3] + 
      	            rnorm(1, mean=0, sd=control_noise)
      state[[tx]]$y[3] <- state[[tx]]$y[3] + 
      	            rnorm(1, mean=0, sd=control_noise)
      tx <- tx + 1
    }
    #cat("\n")

    acc_x <- rep(NA, length(t_vec))
    vel_x <- rep(NA, length(t_vec))
    pos_x <- rep(NA, length(t_vec))
    acc_y <- rep(NA, length(t_vec))
    vel_y <- rep(NA, length(t_vec))
    pos_y <- rep(NA, length(t_vec))
    for( tx in 1:length(t_vec)){ 
      pos_x[tx] <- state[[tx]]$x[1]
      vel_x[tx] <- state[[tx]]$x[2]
      acc_x[tx] <- state[[tx]]$x[3]
      pos_y[tx] <- state[[tx]]$y[1]
      vel_y[tx] <- state[[tx]]$y[2]
      acc_y[tx] <- state[[tx]]$y[3]
    }
  
    traj_mat[n,1,] <- pos_x
    traj_mat[n,2,] <- pos_y
  
    vel_mat[n,1,] <- vel_x
    vel_mat[n,2,] <- vel_y
  }
  return(traj_mat)
} 
