import os
import mmg_coefficients
wave_force = {}
relative_angle = []
def main():
    print("Wave data read code excute")

    now_path = os.getcwd()
    data_path = os.path.join(now_path, 'wave_data')
    velocity = [0.2, 2, 4, 6, 8 ,10, 12]
    global wave_force
    for velocity in velocity:
        data_name = "SWAN_DRIFT_" + str(velocity) + "kts_0.0ms.txt"
        data_file = os.path.join(data_path, data_name)
        data = []
        for line in open(data_file, 'r').readlines():
            data.append(line.split())

        read_line_position = 17
        heading_data = True
        wave_force_velocity = {}

        while read_line_position < len(data):
            if heading_data:
                heading_angle = float(data[read_line_position][3])
                heading_data = False
                wave_force_data = []
                read_line_position += 1
            else:
                X_W = float(data[read_line_position][1])
                Y_W = float(data[read_line_position][2])
                N_W = float(data[read_line_position][-1])
                wave_length = float(data[read_line_position][0])
                wave_force_data_element = [wave_length, X_W, Y_W, N_W]
                wave_force_data.append(wave_force_data_element)

                read_line_position += 1
                if wave_length == 0.850:
                    heading_data = True
                    wave_force_velocity[heading_angle] = list(wave_force_data)
        wave_force[velocity] = wave_force_velocity

    print("Wave data read and make dictionary complete")

    if __name__ == "wave_data_check":
        print("Return wave_force data and this code close")
        return 0
    elif __name__ == "__main__":
        # print(wave_force)
        return 0

def match_wave_force(ship, velocity, ship_heading_angle, wave_angle, wave_length):
    global wave_force
    '''
    wave_force : dict
        dictionary로 1번째 key는 velocity, value는 dictionary
        위의 value dictionary는 key가 heading_angle, value는 list
        list는 wave_length, X_W, Y_W, N_W 순서로 들어가있음
        
        wave_force[velocity][heading_angle][wave_length][X_W, Y_W, N_W] 
    '''
    
    data_velocities = [0.2, 2, 4, 6, 8 ,10, 12]
    data_headings = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
    
    relative_wave_angle = wave_angle + ship_heading_angle
    if relative_wave_angle > 360:
        relative_wave_angle -= 360
    elif relative_wave_angle < 0:
        relative_wave_angle += 360

    # global relative_angle
    # relative_angle.append(relative_wave_angle)
    # print(relative_angle)
    # print(len(relative_angle))
    
    if ship.name == "S175":
        if wave_length == 0.7:
            wave_length = 0.7
        elif wave_length == 1.0:
            wave_length = 0.6
        elif wave_length == 1.2:
            wave_length = 0.55
            
    if velocity in data_velocities:
        if relative_wave_angle in data_headings:
            wave_force_data = wave_force[velocity][relative_wave_angle]
            for i in wave_force_data:
                if i[0] == wave_length:
                    X_W = i[1]
                    Y_W = i[2]
                    N_W = i[3]
                    return X_W, Y_W, N_W
        else:

            lower_heading = max(h for h in data_headings if h < relative_wave_angle)
            upper_heading = min(h for h in data_headings if h > relative_wave_angle)
            
            lower_wave_force_data = wave_force[velocity][lower_heading]
            upper_wave_force_data = wave_force[velocity][upper_heading]


            lower_wave_force = next((i for i in lower_wave_force_data if i[0] == wave_length), None)
            upper_wave_force = next((i for i in upper_wave_force_data if i[0] == wave_length), None)

            if lower_wave_force and upper_wave_force:
                X_W = lower_wave_force[1] + (upper_wave_force[1] - lower_wave_force[1]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                Y_W = lower_wave_force[2] + (upper_wave_force[2] - lower_wave_force[2]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                N_W = lower_wave_force[3] + (upper_wave_force[3] - lower_wave_force[3]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))

                return X_W, Y_W, N_W
            else:
                print("No wave force data available for the given wave length and relative wave angle.")
                return None
    else:
        lower_velocity = max(v for v in data_velocities if v < velocity)
        upper_velocity = min(v for v in data_velocities if v > velocity)

        if relative_wave_angle in data_headings:
            lower_wave_force_data = wave_force[lower_velocity][relative_wave_angle]
            upper_wave_force_data = wave_force[upper_velocity][relative_wave_angle]

            lower_wave_force = next((i for i in lower_wave_force_data if i[0] == wave_length), None)
            upper_wave_force = next((i for i in upper_wave_force_data if i[0] == wave_length), None)

            if lower_wave_force and upper_wave_force:
                X_W = lower_wave_force[1] + (upper_wave_force[1] - lower_wave_force[1]) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))
                Y_W = lower_wave_force[2] + (upper_wave_force[2] - lower_wave_force[2]) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))
                N_W = lower_wave_force[3] + (upper_wave_force[3] - lower_wave_force[3]) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))

                return X_W, Y_W, N_W
            else:
                print("No wave force data available for the given wave length and velocity.")
                return None
        else:
            lower_heading = max(h for h in data_headings if h < relative_wave_angle)
            upper_heading = min(h for h in data_headings if h > relative_wave_angle)

            lower_wave_force_data_lower_heading = wave_force[lower_velocity][lower_heading]
            lower_wave_force_data_upper_heading = wave_force[lower_velocity][upper_heading]
            upper_wave_force_data_lower_heading = wave_force[upper_velocity][lower_heading]
            upper_wave_force_data_upper_heading = wave_force[upper_velocity][upper_heading]

            lower_wave_force_lower_heading = next((i for i in lower_wave_force_data_lower_heading if i[0] == wave_length), None)
            lower_wave_force_upper_heading = next((i for i in lower_wave_force_data_upper_heading if i[0] == wave_length), None)
            upper_wave_force_lower_heading = next((i for i in upper_wave_force_data_lower_heading if i[0] == wave_length), None)
            upper_wave_force_upper_heading = next((i for i in upper_wave_force_data_upper_heading if i[0] == wave_length), None)

            if lower_wave_force_lower_heading and lower_wave_force_upper_heading and upper_wave_force_lower_heading and upper_wave_force_upper_heading:
                X_W_lower = lower_wave_force_lower_heading[1] + (lower_wave_force_upper_heading[1] - lower_wave_force_lower_heading[1]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                Y_W_lower = lower_wave_force_lower_heading[2] + (lower_wave_force_upper_heading[2] - lower_wave_force_lower_heading[2]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                N_W_lower = lower_wave_force_lower_heading[3] + (lower_wave_force_upper_heading[3] - lower_wave_force_lower_heading[3]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))

                X_W_upper = upper_wave_force_lower_heading[1] + (upper_wave_force_upper_heading[1] - upper_wave_force_lower_heading[1]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                Y_W_upper = upper_wave_force_lower_heading[2] + (upper_wave_force_upper_heading[2] - upper_wave_force_lower_heading[2]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))
                N_W_upper = upper_wave_force_lower_heading[3] + (upper_wave_force_upper_heading[3] - upper_wave_force_lower_heading[3]) * ((relative_wave_angle - lower_heading) / (upper_heading - lower_heading))

                X_W = X_W_lower + (X_W_upper - X_W_lower) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))
                Y_W = Y_W_lower + (Y_W_upper - Y_W_lower) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))
                N_W = N_W_lower + (N_W_upper - N_W_lower) * ((velocity - lower_velocity) / (upper_velocity - lower_velocity))

                return X_W, Y_W, N_W
            else:
                print("No wave force data available for the given wave length, velocity and relative wave angle.")
                return None        
        
        
def wave_force_plot():
    import matplotlib.pyplot as plt
    import numpy as np
    import mmg_coefficients

    global wave_force
    wave_force_1 = wave_force[12]
    L = 170
    B = 25.4
    wave_length = [(2*np.pi*9.81/L)/(0.350 + 0.050 * i)**2 for i in range(0, 11)]
    print(wave_length)
    heading_angle = 180
    wave_force_data = wave_force_1[heading_angle]
    ship_type = mmg_coefficients.S175
    X_W = []
    Y_W = []
    N_W = []
    wave_amplitude = 1
    print((1025 * 9.81 * B ** 2 * wave_amplitude ** 2 / L))
    for i in wave_force_data:
        X_W.append(-i[1] * 1025/ (1025 * 9.81 * B ** 2 * wave_amplitude ** 2 / L))
        Y_W.append(-i[2] * 1025/ (1025 * 9.81 * B ** 2 * wave_amplitude ** 2 / L))
        N_W.append(-i[3] * 1025/ (1025 * 9.81 * B ** 3 * wave_amplitude ** 2 / L))

    # print(X_W, Y_W, N_W)
    plt.plot(wave_length, X_W, label = "X_W")
    # plt.plot(wave_length, Y_W, label = "Y_W")
    # plt.plot(wave_length, N_W, label = "N_W")
    plt.legend()
    plt.xlabel("wave_length")
    plt.xlim(0,3)
    plt.ylabel("Wave force [N]")
    plt.ylim(-1,18)
    plt.title("Wave force")
    plt.show()

if __name__ == "__main__":
    main()
    L = 170
    B = 25.4
    # X_W, Y_W, N_W = match_wave_force(mmg_coefficients.S175, 4, 226, 164, 0.7)
    # X_W = - X_W * 1000/ (1025 * 9.81 * B ** 2 * 1 ** 2 / L)
    # Y_W = - Y_W * 1000/ (1025 * 9.81 * B ** 2 * 1 ** 2 / L)
    # N_W = - N_W * 1000/ (1025 * 9.81 * B ** 3 * 1 ** 2 / L)
    # print(X_W, Y_W, N_W)
    wave_force_plot()
elif __name__ == "wave_data_check":
    main()