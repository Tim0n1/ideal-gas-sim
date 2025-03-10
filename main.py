from adodbapi.apibase import pythonTimeConverter
from manim import *
import numpy as np
import time


def elastic_collision(m1, m2, velocity1, velocity2, pos1, pos2):
    r = pos1 - pos2
    n = r / np.linalg.norm(r)
    t = np.array([-n[1], n[0], 0])
    v1n = np.dot(velocity1, n)
    v2n = np.dot(velocity2, n)
    v1t = np.dot(velocity1, t)
    v2t = np.dot(velocity2, t)
    v1n_prime = (m1 - m2) * v1n / (m1 + m2) + 2 * m2 * v2n / (m1 + m2)
    v2n_prime = (m2 - m1) * v2n / (m1 + m2) + 2 * m1 * v1n / (m1 + m2)

    v1t_prime, v2t_prime = v1t, v2t


    v1_prime = v1n_prime * n + v1t_prime * t
    v2_prime = v2n_prime * n + v2t_prime * t

    return v1_prime, v2_prime

def find_grid_cell(x,y, cell_width, start):
    c = int((x - start[0]) // cell_width[0])
    r = int((y - start[1]) // cell_width[1])
    return r,c



class PP(Scene):
    n = 15
    m = 1  # particle mass
    time = 0


    def construct(self):
        left_border = -config.frame_width / 2
        right_border = config.frame_width / 2
        top_border = config.frame_height / 2
        bottom_border = -config.frame_height / 2

        cell_width_x = abs(left_border - right_border) / 16
        cell_width_y = abs(bottom_border - top_border) / 16
        grid = np.zeros(shape=(16,16))

        # Timer
        timer = DecimalNumber(0, num_decimal_places=2)
        timer.to_corner(UP)

        self.add(timer)

        # Draw boundary lines
        line1 = Line(start=(0, top_border, 0), end=(0, bottom_border, 0), stroke_width=5)  # center line
        line2 = Line(start=(left_border, bottom_border, 0), end=(left_border, top_border, 0))
        line3 = Line(start=(right_border, bottom_border, 0), end=(right_border, top_border, 0))
        label_temp1 = Text("Temp: 0.00", font_size=36).to_corner(DOWN + RIGHT)
        label_temp2 = Text("Temp: 0.00", font_size=36).to_corner(DOWN + LEFT)
        entropy_label = Text("Entropy: 0.00", font_size=36).to_corner(UP + RIGHT)
        self.add(line2, line3)

        # Create particles and velocities
        gas1 = VGroup(*[Circle(0.1, color=WHITE, fill_opacity=1).move_to((i, 0, 0)) for i in np.linspace(-6, -0.5, self.n)])
        gas2 = VGroup(*[Circle(0.1, color=RED, fill_opacity=1).move_to((i, 0, 0)) for i in np.linspace(0.5, 6, self.n)])
        velocities = np.random.uniform(-1, 1, size=(2,self.n, 3))
        vel2 = np.random.uniform(-3,3,size=(self.n,3))
        velocities[1] = vel2
        velocities[:,:, 2] = 0

        def update_timer(mob, dt):
            self.time += dt
            mob.set_value(round(self.time,3))

        def particle_motion(group, dt, type=0):
            #print(dt)
            v = velocities[type,:,:]
            for i in range(len(group)):
                group[i].shift(dt * v[i])

        def check_collision(group, type=0):
            collides = {i: False for i in range(len(group))}
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    # Calculate distance between centers of the two particles
                    distance = np.linalg.norm(group[i].get_center() - group[j].get_center())
                    # Check if particles are colliding (based on their radius)
                    if distance < (group[i].get_width() / 2 + group[j].get_width() / 2) + 0.06:
                        group[i].set_color(YELLOW)
                        group[j].set_color(YELLOW)
                        collides[i] = True
                        collides[j] = True
                        u1, u2 = elastic_collision(1, 1, velocities[type,i,:], velocities[type,j,:], group[i].get_center(),
                                                    group[j].get_center())
                        velocities[type,i,:] = u1
                        velocities[type,j,:] = u2

            for i in range(len(group)):
                if not collides[i]:
                    if type == 1:
                        group[i].set_color(RED)
                    else:
                        group[i].set_color(WHITE)

                x, y, _ = group[i].get_center()
                if x - group[i].get_width() <= left_border or x + group[i].get_width() >= right_border:
                    velocities[type,i, 0] *= -1
                if y - group[i].get_width() <= bottom_border or y + group[i].get_width() >= top_border:
                    velocities[type,i, 1] *= -1
                if self.time <= 30:
                    if type == 0:
                        if x + group[i].get_width() >= 0:
                            velocities[type,i, 0] *= -1
                    if type == 1:
                        if x - group[i].get_width() <= 0:
                            velocities[type, i, 0] *= -1

        def label_temp1_updater(label1):
            v_2_mean = (np.linalg.norm(velocities[0], axis=1)**2).mean()
            e_k_mean = v_2_mean / 2
            new_label = Text(f"Temp: {e_k_mean:.2f}", font_size=36).to_corner(LEFT + DOWN)
            label1.become(new_label)

        def label_temp2_updater(label1):
            v_2_mean = (np.linalg.norm(velocities[1], axis=1) ** 2).mean()
            e_k_mean = v_2_mean / 2
            new_label = Text(f"Temp: {e_k_mean:.2f}", font_size=36).to_corner(RIGHT + DOWN)
            label1.become(new_label)

        def line1_updater(mob,dt):
            if self.time > 30:
                mob.set_opacity(0)


        def entropy_updater(e_label, g1, g2):
            grid.fill(0)
            for particle in g1:
                r,c = find_grid_cell(particle.get_x(), particle.get_y(), (cell_width_x, cell_width_y), (left_border, bottom_border))
                grid[r,c] += 1
            for particle in g2:
                r,c = find_grid_cell(particle.get_x(), particle.get_y(), (cell_width_x, cell_width_y), (left_border, bottom_border))
                grid[r,c] += 1
            p_i = grid / (len(g1)+ len(g2))
            p_i = p_i[p_i > 0]
            epsilon = 1e-10
            entropy = -(p_i*np.log(p_i + epsilon)).sum()
            new_label = Text(f"Entropy: {entropy:.5f}", font_size=36).to_corner(UP + RIGHT)
            e_label.become(new_label)

        gas1.add_updater(lambda x, dt: particle_motion(x, dt,0))
        gas1.add_updater(lambda x: check_collision(x,0))
        gas2.add_updater(lambda x, dt: particle_motion(x, dt,1))
        gas2.add_updater(lambda x: check_collision(x,1))
        line1.add_updater(line1_updater)
        label_temp1.add_updater(label_temp1_updater)
        label_temp2.add_updater(label_temp2_updater)
        entropy_label.add_updater(lambda x: entropy_updater(x, gas1, gas2))
        timer.add_updater(update_timer)


        self.add(gas1,gas2, label_temp1, label_temp2, entropy_label, line1)
        self.wait(60)
