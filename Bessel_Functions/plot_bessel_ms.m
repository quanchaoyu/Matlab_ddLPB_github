function plot_bessel_ms
% plot two figures illustrating i_n(x) and k_n(x) 

figure('Position', [900, 200, 700, 280])
subplot(1,2,1)
x = 0:0.01:6;
for n = 0:6
    in = besseli_ms(n,x);
    plot(x,in,'-','LineWidth',1.5)
    hold on
end
axis([0 5 0 10])
grid on
legend('i_0','i_1','i_2','i_3','i_4','i_5','i_6','Location','Best')
%title('Modified Spherical Bessel Functions of the First Kind')
xlabel('x')
ylabel('i_n(x)')

subplot(1,2,2)
x = 0:0.01:6;
for n = 0:6
    kn = besselk_ms(n,x);
    plot(x,kn,'-','LineWidth',1.5)
    hold on
end
axis([0 5 0 10])
grid on
legend('k_0','k_1','k_2','k_3','k_4','k_5','k_6','Location','Best')
%title('Modified Spherical Bessel Functions of the Second Kind')
xlabel('x')
ylabel('k_n(x)')

end