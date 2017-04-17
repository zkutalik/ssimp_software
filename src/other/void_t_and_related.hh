namespace void_t_and_related {

    template<typename ...T>
    struct voider { using type = void; };
    template<typename ...T>
    using void_t = typename voider<T...>:: type;

    template<typename, template<class...> class Template, typename ...Args>
    struct can_apply_impl
        : public std:: false_type {};
    template<template<class...> class Template, typename ...Args>
    struct can_apply_impl<void_t<Template<Args...>>, Template, Args...>
        : public std:: true_type {
            using res = Template<Args...>;
    };
    template<template<class...> class Template, typename ...Args>
    using can_apply = can_apply_impl<void, Template, Args...>;
} // namespace void_t_and_related
