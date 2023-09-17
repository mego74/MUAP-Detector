import matplotlib.pyplot as plt
import numpy as np
import statistics
import math

T = 20
DIFF_TH = 12.75**4.2

number_color_pairs = {
        0: "red",
        1: "green",
        2: "blue",
        3: "yellow",
        4: "orange",
        5: "purple"
    }


def moving_average(data, window_size):
    weights = np.repeat(1.0, window_size) / window_size
    filtered_data = np.convolve(data, weights, 'valid')
    return filtered_data.tolist()


def synchronize(MUAP_Times, filtered_signal):
    for i in range(0, len(MUAP_Times)):
        x = MUAP_Times[i]
        max_index = filtered_signal.index(max(filtered_signal[x:x+T+5]))
        temp = max_index - T//2
        MUAP_Times[i] = temp


def compute_occurrences(filtered_signal, THRESHOLD):
    time_stamp = []
    i = 0
    while i < len(filtered_signal):
        count = 0
        if filtered_signal[i] > THRESHOLD:
            for k in range(i, i+T):
                if filtered_signal[k] > THRESHOLD:
                    count += 1
            if count == 20:
                time_stamp.append(i)
            i = i+T
            j = i
            b = False
            while not b and j < len(filtered_signal)-1:
                if filtered_signal[j+1] > THRESHOLD:
                    i = i+T
                    j = i
                else:
                    b = True
        i = i + 1
    return time_stamp


def compute_summation(M, K):
    # Initialize the summation variable
    D = 0
    for i in range(T):
        D += (M[i] - K[i]) ** 2
    return D


def compute_average(list1, list2):
    averages = []
    for i in range(len(list1)):
        average = (list1[i] + list2[i]) / 2
        averages.insert(i, average)
    return averages


def compute_templates(filtered_signal, first_occurrences_synch, MUAPS_Templates,template_indices):
    templates = []
    start = first_occurrences_synch[0]
    end = first_occurrences_synch[0] + T
    inner_list = filtered_signal[start:end]
    templates.insert(0, inner_list)
    for i in first_occurrences_synch:
        start_index = i
        end_index = i + T
        M = filtered_signal[start_index:end_index]
        detected_template = False
        temp_list = []
        for j in range(0, len(templates)):
            K = templates[j]
            D = compute_summation(M, K)
            temp_list.append(D)
        min_diff = min(temp_list)
        min_diff_index = temp_list.index(min_diff)
        if min_diff < DIFF_TH:
            detected_template = True
            templates[min_diff_index] = compute_average(templates[min_diff_index], M)
            template_indices[min_diff_index].append(start_index)
            MUAPS_Templates[start_index] = min_diff_index
        if not detected_template:
            templates.append(M)
            MUAPS_Templates[start_index] = len(templates)-1
    return templates


def compute_spectrum(signal_size, MUAPS_Templates, template_number):
    MU_Binary_Vectors = [0] * signal_size
    for key in MUAPS_Templates:
        if MUAPS_Templates[key] == template_number:
            MU_Binary_Vectors[key] = 1
    return MU_Binary_Vectors


def main():
    file_path = '/Users/amna_elsaqa/Downloads/Project-3/Data.txt'  # Replace 'filename.txt' with the actual path of your text file
    with open(file_path, 'r') as file:
        signal = file.readlines()
    required_samples = signal
    original_signal = [float(x) for x in required_samples]
    rectified_signal = [abs(x) for x in original_signal]
    filtered_signal = moving_average(rectified_signal, T)
    std_deviation = math.ceil((statistics.stdev(filtered_signal[0:292]) * 3) * 10) / 10
    template_indices = [[], [], [], [], []]
    first_occurrences = compute_occurrences(filtered_signal, std_deviation)
    synchronize(first_occurrences, original_signal)
    first_occurrences.remove(23042)
    first_occurrences.remove(27977)
    first_occurrences.remove(41153)
    MUAPS_Templates = {}
    templates = compute_templates(original_signal, first_occurrences, MUAPS_Templates, template_indices)
    template_indices[1] = template_indices[1] + template_indices[4]
    template_indices[2] = template_indices[2] + template_indices[3]
    template_indices.pop(3)
    template_indices.pop(3)
    temp_1 = compute_average(templates[1], templates[4])
    temp_2 = compute_average(templates[2], templates[3])
    templates[1] = temp_1
    templates[2] = temp_2
    templates.pop(3)
    templates.pop(3)

    for key in MUAPS_Templates:
        if MUAPS_Templates[key] == 4:
            MUAPS_Templates[key] = 1
        if MUAPS_Templates[key] == 3:
            MUAPS_Templates[key] = 2

    x = list(range(len(original_signal[30000:35000])))
    plt.plot(x, original_signal[30000:35000])
    y_line = std_deviation
    plt.plot([min(x), max(x)], [y_line, y_line], 'r--', label='Line at y={}'.format(y_line))
    star_coordinates = []
    for key in MUAPS_Templates:
        if 30000 <= key <= 35000:
            key = key - 30000
            star_coordinates.append((key, 850))
    print(MUAPS_Templates)
    binary_vector1 = compute_spectrum(len(signal), MUAPS_Templates, 0)
    binary_vector2 = compute_spectrum(len(signal), MUAPS_Templates, 1)
    binary_vector3 = compute_spectrum(len(signal), MUAPS_Templates, 2)
    binary_vector = [[], [], []]
    binary_vector[0] = binary_vector1
    binary_vector[1] = binary_vector2
    binary_vector[2] = binary_vector3
    print(binary_vector[0])
    for i in range(len(binary_vector)):
        fft_result = np.fft.fft(binary_vector[i])
        abs_result = np.abs(fft_result)
        binary_vector[i] = abs_result

    for x, y in star_coordinates:
        plt.scatter(x, y, marker='*', s=70, color=number_color_pairs[MUAPS_Templates[x+30000]])
    plt.xlabel('n')
    plt.ylabel('x(n)')
    plt.title('Plot of Data')
    plt.grid(True)
    plt.show()
    for i in range(0,len(templates)):
        x = list(range(len(templates[i])))
        plt.plot(x, templates[i])
        plt.xlabel('n')
        plt.ylabel('x(n)')
        plt.title('Plot of Data')
        plt.grid(True)
        plt.show()
    for i in range(0, len(binary_vector)):
        x = list(range(len(binary_vector[i])))
        plt.plot(x, binary_vector[i])
        plt.xlabel('n')
        plt.ylabel('x(n)')
        plt.title('Plot of Data')
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    main()

